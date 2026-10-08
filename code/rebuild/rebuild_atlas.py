#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
rebuild_atlas.py -- regenerate every derived artefact of the phylogenetic atlas
from the four ground-truth staging files plus the shipped trees.

===============================================================================
WHAT THIS SCRIPT IS
===============================================================================
The 2026-08-20 audit (Paper/audit-2026-08-20/report.pdf; 77 confirmed findings,
23 blocking) established that the atlas's derived ledgers had drifted apart from
each other and from the deposited data, and that the published 23.7 % coverage
headline was a hard-coded literal reachable from no numerator/denominator pair
in the deposit.

The repair discipline is:

    FOUR TRUTH FILES + THE SHIPPED TREES  ->  ONE SCRIPT  ->  EVERY ARTEFACT

    ABSOLUTE RULE: NO REPORTED QUANTITY IS HARD-CODED.  Every number written by
    this script is computed here, at run time, from the truth files and the
    deposited Newick files.  A numeric literal below is a bug unless it is
    (a) a declared RULE PARAMETER (tolerance, size cap, elasticity step),
    (b) a formatting width, or (c) an array index.  All rule parameters live in
    the RULE PARAMETERS block and are echoed into headline_numbers.json.

===============================================================================
THE ADOPTED INCLUSION-CRITERION FRAMEWORK  (the "coverage ladder over one
basis table")
===============================================================================
SCOPE.  Numerator universe = the EUKARYOTIC CANONICAL trees (the canonical
trees of truth_dated.csv, minus the prokaryotic partitions).  Denominator =
described species summed over the EUKARYOTIC partitions of truth_species.csv.
Prokaryotes are reported on their own line and never enter a eukaryotic
numerator or denominator.  The all-life variant is emitted, labelled, for the
SI and for Figure 1 only.

THE ATOM.  Coverage is a predicate on a TIP, never on a partition.  No
partition contributes zero because of what it is; it contributes zero only when
zero of its tips satisfy the predicate.  NO PARTITION SLUG APPEARS ANYWHERE IN
THIS FILE.  Every per-partition difference is read out of the truth files or
measured from the deposited trees.

TWO MEASURED AXES.

  A. Molecular support of placement, f_mol in [0, 1].  For each canonical tree,
     the fraction of its distinct standardised labels that also occur in the
     molecular analysis deposited as that tree's source:

         f_mol = |labels(shipped tree) INTERSECT labels(molecular source)|
                 / |labels(shipped tree)|

     ONE FORMULA FOR EVERYONE.  A tree with no imputation and no grafting has
     the tree itself as its molecular source and scores 1 by the same formula,
     not by exception.  Where no source tree is deposited but the source
     publication states a fraction, f_mol is that fraction, flagged.  Where both
     exist, the MINIMUM is taken (conservative and monotone).

  B. Time calibration, measured from GEOMETRY ONLY.  is_chronogram is the
     root-to-tip-CV verdict of truth_dated.csv, which the audit measured with
     the recursive Newick parser in Paper/audit-2026-08-20/code/common.py.  This
     script IMPORTS that parser (never a reimplementation: a naive one-pass
     parser gets ultrametricity wrong) and RE-VERIFIES every canonical tree's CV
     on each run.  The ledger's `dated` column is never read.  roottotip_cv
     ships alongside the bit so a reader may set their own tolerance.

  C. Uncertainty, interval_recoverable, is archival_class == 'hpd_recoverable'
     in truth_provenance.csv, which is already scoped to the canonical tree.
     divergence_ci_available is deliberately NOT used: it is recorded at the
     wrong grain (it is true for partitions whose canonical tree is not a
     chronogram at all).

THE RUNGS (nested conjuncts, monotone by construction).

    C0  NAMED      tip label occurs in a retained eukaryotic canonical tree
    C1  MOLECULAR  C0 and the tip is in the molecular analysis behind that tree
    C1d DATED      C1 and is_chronogram
    C1b BOUNDED    C1d and interval_recoverable

    Numerator(rung) = SUM over partitions of c0_units * w(rung),
    w in {1, f_mol} x {1, 1[chronogram]} x {1, 1[interval]}.

    There is NO CAP.  No eukaryotic canonical partition has more distinct tip
    labels than its partition has described species; the invariant check proves
    it on every run, and a min(x, 100) cap is therefore not merely unnecessary
    but a mask over a defect if one ever appears.

NAMING DISCIPLINE (a correction to the framework, not a stylistic choice).
C1 is "molecularly supported placement", operationalised by TIP EVIDENCE within
the analysis that produced the tree.  It is NOT a predicate on "placement
traceability": f_mol is an axis-(a) quantity.  The genuine axis-(b) distinction
is reported as a NOMINAL ANNOTATION (topology_assembly) and never as a rung,
because the deposit carries no measured per-tip placement-warrant flag.

MONOTONICITY.  Stated and checked as: the rungs are nested weighted sums with
weights in [0, 1], hence non-increasing; verified per partition on every run.
Nested tip SETS are NOT claimed -- the computed object is an expectation.  Until
every partition is enumerated in molecular_tips.csv, C1 and below are counted in
TIP-EQUIVALENTS, never in species; headline_numbers.json labels them so.

DATING IS A FILTER, NEVER A COVERAGE COLUMN in the per-partition table.  The
per-partition table carries C0 %, C1 %, is_chronogram, roottotip_cv,
interval_kind.  A completely sampled but uncalibrated partition therefore reads
"91.39 / 91.39 / not ultrametric (CV 0.296) / --", which is the true statement,
instead of "91.39 / 91.39 / 0 / 0", which is not.  The nested chain survives in
the abstract sentence, where it is a refinement of one proposition.

TOPOLOGY ASSEMBLY is a NOMINAL annotation over
{primary_analysis, compiled, automated_synthesis}, a pure function of two
deposited fields of the PRE-AUDIT ledger Data/tree_provenance.csv.  It has no
ordering, no ">=" operator and no strong/weak reading, and Methods must say that
it is derived from UNRECONCILED ledger fields and is therefore an annotation,
not a criterion.  It is reported two ways: the Figure 2 scatter, and one
leave-class-out sensitivity line.

===============================================================================
INPUTS
===============================================================================
REQUIRED -- the script REFUSES TO RUN if any is missing:
  Data/_rebuild/truth_dated.csv        one row per shipped .nwk; the geometry
                                       columns are re-measured on every run
  Data/_rebuild/truth_tips.csv         one row per shipped .nwk; likewise
  Data/_rebuild/truth_species.csv       62 rows, one per partition
  Data/_rebuild/truth_provenance.csv    26 canonical + alternative rows
  Data/standardized/dictionary.csv     numeric tip id -> standardised name
  Data/standardized/full_mapping.csv   original label -> standardised name
  Data/standardized/trees/*.nwk        the shipped trees
  Paper/audit-2026-08-20/code/common.py   the reference recursive Newick parser

READ-ONLY, and demoted in Methods to "the pre-audit ledger":
  Data/tree_provenance.csv   ONLY `inference`, `otl_pipeline_involvement` and
                             `tree_class` are read, and only to populate the
                             nominal topology_assembly annotation and the S5
                             tree-class column.

PASSTHROUGH (fields the truth files deliberately do not carry; each is tagged
PASSTHROUGH in the run report so a reader can see it was not recomputed):
  estimated_total, confidence   curated info.yaml -> data_estimates.csv
                                -> existing partition_summary.csv (majority vote)
  S7 literature rows            the published-estimate compilation inlined in
                                Paper/shared/supplementary_body.tex.  Only the
                                "This study" row of S7 is recomputed.

===============================================================================
OUTPUTS  (all written atomically: temp file in the same directory, then rename)
===============================================================================
THE DEPOSIT CONTRACT
  Data/coverage_basis.csv        62 rows, ONE PER PARTITION, including the 36
                                 with no tree (labels_unique = 0) and the
                                 prokaryotic ones (flagged).  Every numerator
                                 and every denominator in the paper is a column
                                 sum of THIS file, which structurally forecloses
                                 a phantom denominator forever.
  Data/molecular_tips.csv        (partition, standardized_id) for every
                                 partition whose molecular source is a deposited
                                 tree -- so C1 is a genuine enumerable,
                                 dedupable, prunable tip set there.
  Data/atlas_coverage.py         the ~30-line reader.  Every published figure is
  Data/atlas_coverage.R          one call of it.  Shipped in Data Availability.

DERIVED LEDGERS AND FIGURE INPUTS
  Data/site/partition_summary.csv          62 rows, web-app schema preserved
  Paper/inputs/coverage_per_partition.csv  the per-partition ladder table
  Paper/inputs/fig1_coverage_bars.csv      Figure 1 (all-life inventory)
  Paper/inputs/fig2_coverage_scatter.csv   Figure 2 (f_mol x C0 x assembly x size)
  Paper/inputs/fig3_c0c1_bars.csv          Figure 3 (C0/C1 bars + structural zeros)
  Paper/inputs/fig3_dated_undated_counts.csv
  Paper/inputs/table_s5_tree_provenance.csv / .tex
  Paper/inputs/table_s6_figure1_source_data.csv / .tex
  Paper/inputs/table_s6_alternative_trees.csv
  Paper/inputs/table_s7_richness_estimates.csv / .tex
  Paper/inputs/headline_numbers.json       THE single source both the figures and
                                           the manuscript read.  Contains every
                                           rung, every sensitivity, every
                                           disclosure and every open flag.

RETIRED, and never emitted by this script:
  the 23.7 % coverage headline (hard-coded at
  Paper/code/figures/fig2_sensitivity_bounds.py:61; implied denominator
  2,084,409 exists in no file) -- WITHDRAWN, not adjusted; and the old middle
  rung "C1 = 8.65 %", which never had a principled definition.  This script
  does not own Paper/code/; the deletion of that file and of its literal is
  reported, not performed.

NOT REBUILT (source manifests and literature inventories, not quantities
derivable from the truth files) -- listed at the end of every run.

===============================================================================
SAFETY PROPERTIES
===============================================================================
  * REFUSES TO RUN if any required input is missing or has zero data rows.
  * ARCHIVES FIRST: every file about to be overwritten is copied into
    Data/_archive/2026-08-20-pre-audit/<same relative path> if and only if no
    copy is already there.  An archived copy is NEVER overwritten -- the
    pre-audit snapshot must keep its deliberately-bad values.
  * ATOMIC WRITES via <name>.tmp-<pid> + os.replace(), so a crash or a full disk
    cannot leave a truncated ledger behind.
  * IDEMPOTENT: byte-identical output on a second run.  Nothing is timestamped,
    nothing is randomised, every collection is sorted by an explicit key, and
    floats go through one fixed formatter.
  * CONCURRENCY: other agents work in this repository simultaneously.  This
    script writes ONLY the files in OUT.  The truth files, Data/tree_provenance.csv,
    Paper/shared/*.tex, Paper/code/ and the pre-audit archive are opened
    READ-ONLY.  Anything else that needs changing is REPORTED, never touched.

USAGE
    python3 scripts/rebuild_atlas.py                # rebuild everything
    python3 scripts/rebuild_atlas.py --dry-run      # compute + diff, write nothing
    python3 scripts/rebuild_atlas.py --only headline_numbers.json
    python3 scripts/rebuild_atlas.py --quiet        # summary lines only
    python3 scripts/rebuild_atlas.py --no-verify-trees   # skip CV re-measurement

Exit codes: 0 success; 1 refusal (missing input, failed invariant); 2 write failure.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import difflib
import importlib.util
import io
import json
import math
import os
import re
import shutil
import sys
from collections import Counter, OrderedDict, defaultdict

# CSV cells in this deposit include very long curatorial notes and, in the
# standardised layer, very long label strings.
csv.field_size_limit(1 << 30)

# The reference Newick parser is recursive; deep ladder-like clades in the large
# trees can exceed the default limit.  Raised deliberately, with a clear error.
sys.setrecursionlimit(300000)


# ==========================================================================
# PATHS.  Resolved relative to the repository root (the parent of scripts/).
# No absolute path is baked in.
# ==========================================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir))

ARCHIVE_DIR = os.path.join(ROOT, "Data", "_archive", "2026-08-20-pre-audit")
REBUILD_DIR = os.path.join(ROOT, "Data", "_rebuild")

TRUTH_FILES = OrderedDict([
    ("dated",      os.path.join(REBUILD_DIR, "truth_dated.csv")),
    ("tips",       os.path.join(REBUILD_DIR, "truth_tips.csv")),
    ("species",    os.path.join(REBUILD_DIR, "truth_species.csv")),
    ("provenance", os.path.join(REBUILD_DIR, "truth_provenance.csv")),
])

STD_DIR = os.path.join(ROOT, "Data", "standardized")
DICTIONARY = os.path.join(STD_DIR, "dictionary.csv")
FULL_MAPPING = os.path.join(STD_DIR, "full_mapping.csv")
STD_TREES = os.path.join(STD_DIR, "trees")
RAW_DIR = os.path.join(ROOT, "Data", "raw")
# Per-tip molecular evidence for partitions whose molecular analysis is not
# deposited as a tree with the shipped labels (basis per_tip_evidence).
MOL_TIP_EVIDENCE = os.path.join(ROOT, "Data", "molecular_tip_evidence.csv")
MOL_EVIDENCE_KINDS = ("source_tree_label", "reference_genome",
                      "analysis_terminal")
# S05 / D18: per (partition, tip) resolution against the Catalogue of Life
# base release, written by scripts/match_col.py.  C0 counts, per eukaryotic
# partition, the distinct accepted names its tips resolve to.
COL_RESOLUTION = os.path.join(STD_DIR, "col_resolution.csv")
COL_RESOLUTION_RESOLVED = ("accepted", "provisional", "synonym",
                           "accepted_orthographic",
                           "provisional_orthographic",
                           "synonym_orthographic")
# S09: the same resolution over every non-canonical eukaryotic tree.
COL_RESOLUTION_NONCANON = os.path.join(STD_DIR,
                                       "col_resolution_noncanonical.csv")
# D14: which ultrametric canonical files are in absolute time.
TIME_SCALE_CHECKS = os.path.join(ROOT, "Data", "time_scale_checks.csv")
_TIME_SCALES = None


def time_scales():
    global _TIME_SCALES
    if _TIME_SCALES is None:
        _TIME_SCALES = {}
        if os.path.exists(TIME_SCALE_CHECKS):
            with open(TIME_SCALE_CHECKS, newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh):
                    _TIME_SCALES[r["filename"].strip()] = \
                        r["time_scale"].strip()
    return _TIME_SCALES


# S02: crown-age intervals measured from the deposits by
# scripts/crown_intervals.py.  The crown POINT is the root height of the
# canonical file, measured here; the interval, and whether one is recoverable
# at all (rung C1b), come from this register and never from a ledger column.
CROWN_INTERVAL_SOURCES = os.path.join(ROOT, "Data", "crown_interval_sources.csv")
CROWN_INTERVALS = os.path.join(ROOT, "Data", "crown_intervals.csv")
_CROWN_INTERVALS = None
# |focal-node height in a summary-tree deposit - canonical root height| above
# which the focal node is not the canonical crown (Myr).
CROWN_POINT_TOL = 1e-3
CROWN_DP = 4        # decimal places for ages in Ma


def crown_intervals():
    """{partition: row} of Data/crown_intervals.csv (all kinds)."""
    global _CROWN_INTERVALS
    if _CROWN_INTERVALS is None:
        _CROWN_INTERVALS = {}
        if os.path.exists(CROWN_INTERVALS):
            with open(CROWN_INTERVALS, newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh):
                    if r["partition"] in _CROWN_INTERVALS:
                        die("%s has two rows for %s"
                            % (rel(CROWN_INTERVALS), r["partition"]))
                    _CROWN_INTERVALS[r["partition"]] = r
    return _CROWN_INTERVALS


# Root height of every shipped file, filled by measure_file().
ROOT_HEIGHT = {}


# Revision 2026-10 check registers (decisions D7 and D8).
OFFTARGET_CHECKS = os.path.join(ROOT, "Data", "offtarget_checks.csv")
PRUNED_TIPS = os.path.join(ROOT, "Data", "pruned_tips.csv")
DATING_DEPOSIT_CHECKS = os.path.join(ROOT, "Data", "dating_deposit_checks.csv")

AUDIT_COMMON = os.path.join(ROOT, "Paper", "audit-2026-08-20", "code", "common.py")

sys.path.insert(0, SCRIPT_DIR)
from canonical_trees import species_key as SPECIES_KEY  # noqa: E402  (D4)

# Read-only passthrough / pre-audit ledger sources.
TREE_PROVENANCE = os.path.join(ROOT, "Data", "tree_provenance.csv")
DATA_ESTIMATES = os.path.join(ROOT, "Data", "data_estimates.csv")
OLD_PARTITION_SUMMARY = os.path.join(ROOT, "Data", "site", "partition_summary.csv")
PARTITION_INFO_DIR = os.path.join(ROOT, "Data", "site", "data", "partitions")
# The pre-revision supplement, archived: the literature rows of the
# richness-estimate table are read from it (build_table_s7_*).
SUPPLEMENTARY_BODY = os.path.join(ROOT, "Paper", "archived", "2026-08-manuscript",
                                  "supplementary_body.tex")

OUT = OrderedDict([
    # ---- the deposit contract ----
    ("coverage_basis.csv",     os.path.join(ROOT, "Data", "coverage_basis.csv")),
    ("molecular_tips.csv",     os.path.join(ROOT, "Data", "molecular_tips.csv")),
    ("atlas_coverage.py",      os.path.join(ROOT, "Data", "atlas_coverage.py")),
    ("atlas_coverage.R",       os.path.join(ROOT, "Data", "atlas_coverage.R")),
    # ---- derived ledgers and figure inputs ----
    ("partition_summary.csv",
     os.path.join(ROOT, "Data", "site", "partition_summary.csv")),
    ("coverage_per_partition.csv",
     os.path.join(ROOT, "Paper", "inputs", "coverage_per_partition.csv")),
    ("fig1_coverage_bars.csv",
     os.path.join(ROOT, "Paper", "inputs", "fig1_coverage_bars.csv")),
    ("fig2_coverage_scatter.csv",
     os.path.join(ROOT, "Paper", "inputs", "fig2_coverage_scatter.csv")),
    ("fig3_c0c1_bars.csv",
     os.path.join(ROOT, "Paper", "inputs", "fig3_c0c1_bars.csv")),
    ("fig3_dated_undated_counts.csv",
     os.path.join(ROOT, "Paper", "inputs", "fig3_dated_undated_counts.csv")),
    ("table_s5_tree_provenance.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_s5_tree_provenance.csv")),
    ("table_s5_tree_provenance.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_s5_tree_provenance.tex")),
    ("table_s6_figure1_source_data.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_s6_figure1_source_data.csv")),
    ("table_s6_figure1_source_data.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_s6_figure1_source_data.tex")),
    ("table_s6_alternative_trees.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_s6_alternative_trees.csv")),
    ("table_s7_richness_estimates.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_s7_richness_estimates.csv")),
    ("table_s7_richness_estimates.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_s7_richness_estimates.tex")),
    ("c0_name_resolution.csv",
     os.path.join(ROOT, "Paper", "inputs", "c0_name_resolution.csv")),
    ("crown_ages.csv",
     os.path.join(ROOT, "Paper", "inputs", "crown_ages.csv")),
    ("opentree_comparison.csv",
     os.path.join(ROOT, "Paper", "inputs", "opentree_comparison.csv")),
    ("table_datasets.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_datasets.csv")),
    ("table_datasets.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_datasets.tex")),
    ("table_datasets_macros.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_datasets_macros.tex")),
    ("crown_ages_macros.tex",
     os.path.join(ROOT, "Paper", "inputs", "crown_ages_macros.tex")),
    ("numbers_macros.tex",
     os.path.join(ROOT, "Paper", "inputs", "numbers_macros.tex")),
    ("table_which_tree.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_which_tree.csv")),
    ("table_which_tree.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_which_tree.tex")),
    ("table_rungs_sensitivity.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_rungs_sensitivity.csv")),
    ("table_rungs_sensitivity.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_rungs_sensitivity.tex")),
    ("table_worked_example.csv",
     os.path.join(ROOT, "Paper", "inputs", "table_worked_example.csv")),
    ("table_worked_example.tex",
     os.path.join(ROOT, "Paper", "inputs", "table_worked_example.tex")),
    ("fig1_coverage_cascade.csv",
     os.path.join(ROOT, "Paper", "inputs", "fig1_coverage_cascade.csv")),
    ("headline_numbers.json",
     os.path.join(ROOT, "Paper", "inputs", "headline_numbers.json")),
])

NOT_REBUILT = OrderedDict([
    ("Paper/inputs/fig2_tree_manifest.csv",
     "source-file manifest; lists raw deposit paths, not a derived quantity"),
    ("Paper/inputs/fig_crown_age_points.csv",
     "written by Paper/code/figures/fig2_crown_age.R from crown_ages.csv; "
     "verify_atlas.py (CROWNREG) compares it with crown_ages.csv"),
    ("Paper/inputs/fig2_species_accumulation / fig4_*.csv",
     "geological time-scale constants and figure-specific source manifests"),
    ("Paper/inputs/figN_atlas_growth.csv",
     "canonical-succession time series; needs a canonical-succession ledger "
     "that is not among the truth files"),
    ("Paper/inputs/data_provenance.csv",
     "per-source provenance ledger; literature metadata, not a derived quantity"),
    ("Paper/inputs/vertebrates_plants_cr_manifest.csv",
     "figure-specific source manifest"),
])


# ==========================================================================
# RULE PARAMETERS.  These are RULES, not reported quantities.  Every one of
# them is echoed into headline_numbers.json under `rule_parameters` so that no
# consumer has to read this file to know what tolerance produced a number.
# ==========================================================================

YES = {"yes", "true", "1"}

# --- molecular-source resolution -----------------------------------------
# A deposited tree is accepted as a partition's molecular source only if the
# intersection it implies reproduces the numerator the ground truth recorded,
# to within this relative tolerance, AND the source's label set is (near) a
# subset of the shipped tree's.
SOURCE_NUMERATOR_TOL = 0.02      # |measured - recorded| / recorded
SOURCE_NONSUBSET_TOL = 0.05      # |source \ shipped| / |source|
# Directory scans (used only when the ground truth names a raw path that is not
# a usable tree file) skip files above this size: the deposit contains
# half-gigabyte posterior tree sets whose first tree adds nothing.
SOURCE_SCAN_MAX_BYTES = 64 * 1024 * 1024
SOURCE_TREE_EXT = (".tre", ".nwk", ".newick", ".tree", ".taxa.txt")
# A '.taxa.txt' source is the list of taxa with sequence data in a deposited
# alignment (Data/molecular_sources/), one label per line, '#' comments.
TAXA_LIST_EXT = ".taxa.txt"

# --- ultrametricity re-verification --------------------------------------
# truth_dated.csv already carries the audit's verdict.  Each run re-measures the
# canonical trees with the audit's own parser and warns if the CV differs by
# more than this, so ground truth cannot silently rot.
CV_RECHECK_TOL = 1e-3
# The chronogram rule: root-to-tip CV at or below this => chronogram.  The same
# value is declared independently in scripts/verify_atlas.py.  Used to REFRESH
# the measured columns of truth_dated.csv from the shipped files on every run
# (refresh_measured_truth), so that no CV and no is_chronogram verdict is ever
# a typed value.
CHRONOGRAM_TOLERANCE = 0.05

# --- denominator elasticity demonstration --------------------------------
# The single largest eukaryotic described-species count is perturbed by this
# fraction, in both directions, to show how much of every printed decimal rests
# on one round number.  A DEMONSTRATION OF ELASTICITY, NOT A RICHNESS CLAIM.
ELASTICITY_FRACTION = 0.20

# --- display / formatting -------------------------------------------------
PCT_DP = 4          # decimal places for percentages in CSVs
FRAC_DP = 6         # decimal places for fractions (f_mol, CV)
# Rung numerators are WEIGHTED sums and are genuinely fractional.  They are
# rounded once, to this many places, and every percentage is formed from the
# rounded value -- so every printed percentage is reachable from a numerator
# and a denominator printed beside it.  Publishing only the integer would
# leave percentages that no printed pair reproduces, which is the defect the
# retired 23.7 % headline was.
NUMERATOR_DP = 4
CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods",
                  "Other animals", "Microbes & protists"]

# --- the nominal topology-assembly classifier (Methods prints this verbatim) --
# A pure function of two fields of the PRE-AUDIT ledger Data/tree_provenance.csv.
# NOMINAL: no ordering, no ">=", no strong/weak reading.
ASSEMBLY_AUTOMATED = "automated_synthesis"
ASSEMBLY_COMPILED = "compiled"
ASSEMBLY_PRIMARY = "primary_analysis"
ASSEMBLY_OTL_TRIGGER = "full_pipeline"
ASSEMBLY_INFERENCE_TRIGGERS = ("synthesis", "literature_compilation",
                               "supertree_re_dated")


def classify_assembly(otl, inference):
    """topology_assembly, verbatim as Methods states it."""
    if (otl or "").strip() == ASSEMBLY_OTL_TRIGGER:
        return ASSEMBLY_AUTOMATED
    if (inference or "").strip() in ASSEMBLY_INFERENCE_TRIGGERS:
        return ASSEMBLY_COMPILED
    return ASSEMBLY_PRIMARY


# --- the rungs ------------------------------------------------------------
RUNGS = OrderedDict([
    ("C0",  "NAMED: the tip label occurs in a retained eukaryotic canonical tree"),
    ("C1",  "MOLECULAR: C0 and the tip is in the molecular analysis behind that tree"),
    ("C1d", "DATED: C1 and the shipped canonical tree is a measured chronogram"),
    ("C1b", "BOUNDED: C1d and the divergence-time interval is recoverable from "
            "the source deposit"),
])


# ==========================================================================
# small utilities
# ==========================================================================

def die(msg, code=1):
    sys.stderr.write("REFUSING TO RUN: " + msg + "\n")
    sys.exit(code)


def truthy(v):
    return str(v).strip().lower() in YES


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
    """Fixed-precision float rendering: the single formatter used everywhere.
    One formatter is an IDEMPOTENCY requirement, not a cosmetic choice."""
    if x is None:
        return ""
    return "{0:.{1}f}".format(float(x), dp)


def rel(path):
    try:
        return os.path.relpath(path, ROOT)
    except ValueError:
        return path


def read_csv(path, required=True):
    if not os.path.exists(path):
        if required:
            die("missing required input: %s" % rel(path))
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s).lower()).strip("_")


def tex_escape(s):
    """Escape a plain string for LaTeX.  Applied to every cell of every table
    emitted here, so a stray '&' or '_' in a citation key cannot break the
    supplement build."""
    if s is None:
        return ""
    out = []
    for ch in str(s):
        if ch == "\\":
            out.append(r"\textbackslash{}")
        elif ch in "&%$#_{}":
            out.append("\\" + ch)
        elif ch == "~":
            out.append(r"\textasciitilde{}")
        elif ch == "^":
            out.append(r"\textasciicircum{}")
        else:
            out.append(ch)
    return "".join(out)


def median(xs):
    """Median of a list of numbers (mean of the two middle values)."""
    s_ = sorted(xs)
    k = len(s_) // 2
    return s_[k] if len(s_) % 2 else (s_[k - 1] + s_[k]) / 2.0


def thousands(n):
    """LaTeX-safe thousands separator matching the supplement's house style."""
    if n is None or n == "":
        return "--"
    try:
        return "{:,}".format(int(n)).replace(",", "{,}")
    except (TypeError, ValueError):
        return tex_escape(n)


def commas(n):
    if n is None or n == "":
        return "--"
    return "{:,}".format(int(round(float(n))))


def pct(num, den):
    """The ONLY place a percentage is formed.  Returns None on a bad
    denominator rather than raising or silently emitting zero."""
    if den in (None, 0):
        return None
    return 100.0 * float(num) / float(den)


# ==========================================================================
# THE REFERENCE NEWICK PARSER -- imported, never reimplemented
# ==========================================================================

def load_audit_parser():
    """Import Paper/audit-2026-08-20/code/common.py and return parse_newick.

    The audit's parser is RECURSIVE-DESCENT and is the one that gets
    ultrametricity right; a naive one-pass parser does not.  The repair brief is
    explicit -- reuse it, do not reinvent it -- so this script refuses to run
    without it rather than falling back on a parser of its own.
    """
    if not os.path.exists(AUDIT_COMMON):
        die("the audit's reference Newick parser is missing:\n"
            "                %s\n"
            "                This script imports parse_newick from it rather "
            "than reimplementing it,\n"
            "                because a naive one-pass parser gets "
            "ultrametricity wrong." % rel(AUDIT_COMMON))
    spec = importlib.util.spec_from_file_location("_atlas_audit_common",
                                                  AUDIT_COMMON)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:                              # noqa: BLE001
        die("could not import the audit parser from %s: %s"
            % (rel(AUDIT_COMMON), exc))
    if not hasattr(mod, "parse_newick"):
        die("%s does not define parse_newick." % rel(AUDIT_COMMON))
    return mod.parse_newick


PARSE_NEWICK = None          # bound in main(), after the refusal checks


# -- WHITESPACE NORMALISATION -------------------------------------------------
# The audit's reference parser is recursive-descent and correct on ultrametry,
# which is why it is imported rather than reimplemented.  It does NOT, however,
# skip whitespace at structural positions: after a comma it tests s[pos] == "("
# directly, so a deposited tree written as ", (" is read as an empty-named leaf
# and the remainder of the file is mis-nested.  Exactly one shipped tree is
# written that way -- cetaceans.nwk, whose 85 leaves the reference parser reads
# as 4 -- and truth_tips.csv's 85 is the correct value.
#
# The fix belongs HERE, not in the audit's frozen common.py: the parser is
# reused unmodified and handed text that has been normalised first.  The
# normalisation is quote-aware and removes whitespace ONLY where it is adjacent
# to a structural character, so a quoted label such as 'Acanthella acuta' -- of
# which the raw deposit has many -- survives intact.  Verified over all 264
# shipped trees: the parse changes for cetaceans.nwk alone, and after
# normalisation every tree agrees with truth_tips.csv exactly.
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
                i = j                      # structural whitespace -> drop
                continue
            out.append(" ")                # whitespace inside a label -> keep
            i = j
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def tree_facts(text):
    """One parse -> (leaf labels, root-to-tip depths).

    Uses the audit's recursive parser, and the audit's own NaN convention:
    a NaN branch length contributes zero to the depth.  Returning both in a
    single traversal is what lets every run re-verify truth_dated.csv's
    ultrametricity verdict without parsing the trees twice.
    """
    try:
        root = PARSE_NEWICK(normalize_newick(text))
    except RecursionError:
        raise RecursionError(
            "the recursive reference parser exhausted the stack; raise "
            "sys.setrecursionlimit further or report the offending tree")
    labels, depths = [], []
    stack = [(root, 0.0)]
    while stack:
        (kids, name, bl), acc = stack.pop()
        if bl != bl:                        # NaN branch length
            bl = 0.0
        d = acc + bl
        if kids:
            for k in kids:
                stack.append((k, d))
        else:
            depths.append(d)
            if name:
                labels.append(name)
    return labels, depths


def coefficient_of_variation(xs):
    if not xs:
        return None
    mean = sum(xs) / len(xs)
    if mean == 0:
        return None
    var = sum((x - mean) ** 2 for x in xs) / len(xs)
    return math.sqrt(var) / abs(mean)


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
# MEASURED TRUTH -- the geometry columns of truth_dated.csv / truth_tips.csv
# ==========================================================================
#
# truth_dated.csv and truth_tips.csv carry two kinds of column.  The MEASURED
# ones (leaves, distinct labels, NaN branch lengths, root-to-tip CV, the
# chronogram verdict, the adopted tip count and its delta to the published
# count) are properties of the shipped files and are re-measured here on every
# run, from the files, with the audit's parser.  The CURATED ones (partition,
# canonical flag, the pre-audit ledger flags, published_ntips and the
# convention text) are carried through untouched.
#
# A shipped tree that has no row is a refusal: its curated columns cannot be
# invented.  A row whose file is no longer shipped is dropped and reported,
# which is how a withdrawn file leaves every count that follows from it.

DATED_MEASURED = ("leaves", "distinct_labels", "nan_branch_lengths",
                  "roottotip_cv", "is_chronogram")
TIPS_MEASURED = ("raw_leaves", "distinct_labels", "tips_adopted",
                 "source_ntips", "delta")


def measure_file(path, record_root=True):
    """(n_leaves, n_distinct, n_nan, cv) on COMMON.leaf_depths' convention:
    every leaf label enters the distinct count, a NaN branch length counts and
    contributes zero depth.  record_root=False for a file that is not a
    shipped tree (a canonical source), whose basename may equal one."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        root = PARSE_NEWICK(normalize_newick(fh.read()))
    depths, labels, nan_bl = [], [], 0
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
            labels.append(name)
    if record_root:
        ROOT_HEIGHT[os.path.basename(path)] = max(depths) if depths else None
    return len(depths), len(set(labels)), nan_bl, coefficient_of_variation(depths)


def refresh_measured_truth(dry_run, warnings):
    """Rewrite the measured columns of truth_dated.csv and truth_tips.csv from
    the shipped files.  Returns {"dated": rows, "tips": rows} for Truth()."""
    on_disk = sorted(fn for fn in os.listdir(STD_TREES) if fn.endswith(".nwk"))
    out = {}
    for key, cols in (("dated", DATED_MEASURED), ("tips", TIPS_MEASURED)):
        path = TRUTH_FILES[key]
        with open(path, newline="", encoding="utf-8") as fh:
            rdr = csv.DictReader(fh)
            header = list(rdr.fieldnames)
            rows = list(rdr)
        out[key] = (path, header, rows)

    have = {r["filename"] for r in out["dated"][2]}
    missing = [fn for fn in on_disk if fn not in have]
    if missing:
        die("shipped tree(s) with no row in truth_dated.csv: %s. Their curated "
            "columns (partition, canonical flag) cannot be invented; add the "
            "rows first." % ", ".join(missing))

    facts = {}
    for fn in on_disk:
        facts[fn] = measure_file(os.path.join(STD_TREES, fn))
    # source_ntips: the leaf count of the source tree a canonical file is built
    # from (Data/canonical_sources.csv), before D7 pruning and D4 collapsing.
    # It is a separate column so that published_ntips keeps the count the
    # published table carried.
    src_leaves = {}
    canon_src = os.path.join(ROOT, "Data", "canonical_sources.csv")
    if os.path.exists(canon_src):
        for r in read_csv(canon_src):
            sp = os.path.join(ROOT, r["source"].strip())
            if not os.path.exists(sp):
                die("Data/canonical_sources.csv names %s as the source of %s, "
                    "but it is not present" % (r["source"], r["filename"]))
            src_leaves[r["filename"].strip()] = measure_file(sp, record_root=False)[0]

    result = {}
    for key, (path, header, rows) in out.items():
        kept = []
        for r in rows:
            fn = r["filename"]
            if fn not in facts:
                warnings.append(
                    "%s: %s is no longer among the shipped trees; its row is "
                    "dropped and every count that follows from the shipped set "
                    "is recomputed without it." % (rel(path), fn))
                continue
            n, nd, nan, cv = facts[fn]
            r = OrderedDict((h, r.get(h, "")) for h in header)
            if key == "dated":
                r["leaves"] = str(n)
                r["distinct_labels"] = str(nd)
                r["nan_branch_lengths"] = str(nan)
                r["roottotip_cv"] = fnum(cv, FRAC_DP) if cv is not None else ""
                r["is_chronogram"] = ("yes" if cv is not None
                                      and cv <= CHRONOGRAM_TOLERANCE else "no")
            else:
                conv = (r.get("convention") or "").split(";")[0].strip()
                if conv == "raw_leaves":
                    adopted = n
                elif conv == "distinct_labels":
                    adopted = nd
                else:
                    die("%s: %s carries unknown tip convention %r"
                        % (rel(path), fn, conv))
                r["raw_leaves"] = str(n)
                r["distinct_labels"] = str(nd)
                r["tips_adopted"] = str(adopted)
                if "source_ntips" in r:
                    r["source_ntips"] = (str(src_leaves[fn]) if fn in src_leaves
                                         else "")
                pub_n = as_int(r.get("published_ntips"))
                r["delta"] = "" if pub_n is None else str(adopted - pub_n)
            kept.append(r)
        text = render_csv(header, kept)
        status = atomic_write(path, text, dry_run, [], warnings)
        if status not in ("unchanged",):
            print("  %-42s %s (measured columns refreshed from the files)"
                  % (rel(path), status))
        result[key] = [dict(r) for r in kept]
    return result


# ==========================================================================
# atomic write + archive-before-overwrite
# ==========================================================================

def archive_original(dest_path, log, warnings=None):
    """Preserve the pre-correction copy of a file we are about to overwrite.

    THE ARCHIVE IS SEALED.  Data/_archive/2026-08-20-pre-audit/ is a ONE-SHOT
    snapshot taken before the corrections were applied; its MANIFEST.md
    enumerates it and names the deliberately-bad values it preserves on
    purpose.  This script must never modify it AND MUST NEVER EXTEND IT.

    Extending it is not harmless: on the second run, an artefact THIS SCRIPT
    produced on the first run would be copied in under the name of a
    "pre-audit" file, silently contaminating the only route back to the
    pre-correction state with a post-correction artefact.  That happened once
    (headline_numbers.json) and is what this guard exists to prevent.

    So: a file that already has an archived copy is left alone, and a file that
    does not is NOT added -- its absence from the snapshot means it was not
    part of the pre-audit state.  If such a file is nonetheless being
    overwritten, that is surfaced as a warning for a human to rule on rather
    than resolved by writing into the sealed archive.
    """
    if not os.path.exists(dest_path):
        return None                       # nothing to preserve; new artefact
    archived = os.path.join(ARCHIVE_DIR, rel(dest_path))
    if os.path.exists(archived):
        return "already-archived"
    if warnings is not None:
        warnings.append(
            "%s is being overwritten but has no copy in the sealed pre-audit "
            "snapshot. It is NOT added to the archive: the snapshot is a "
            "one-shot pre-correction record and adding a post-correction "
            "artefact to it would destroy that guarantee. If this file did "
            "hold pre-audit content, snapshot it by hand before rerunning."
            % rel(dest_path))
    return "not-in-sealed-snapshot"


def atomic_write(dest_path, text, dry_run, log, warnings=None):
    """Write `text` to `dest_path` atomically, archiving any original first."""
    existed = os.path.exists(dest_path)
    if existed:
        with open(dest_path, "r", encoding="utf-8", errors="replace") as fh:
            if fh.read() == text:
                return "unchanged"
    if dry_run:
        return "would-write" if existed else "would-create"

    archive_original(dest_path, log, warnings)
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    tmp = "%s.tmp-%d" % (dest_path, os.getpid())
    try:
        with open(tmp, "w", encoding="utf-8", newline="") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, dest_path)        # atomic on POSIX
    except Exception as exc:              # noqa: BLE001 - re-raised as exit 2
        if os.path.exists(tmp):
            try:
                os.unlink(tmp)
            except OSError:
                pass
        sys.stderr.write("FAILED to write %s: %s\n" % (rel(dest_path), exc))
        sys.exit(2)
    return "written" if existed else "created"


def render_csv(header, rows):
    """Serialise rows to a CSV string with LF line endings.  Deterministic by
    construction: the caller supplies both column order and row order."""
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=header, lineterminator="\n",
                       extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in header})
    return buf.getvalue()


# ==========================================================================
# LOAD + VALIDATE the four truth files
# ==========================================================================

class Truth(object):
    """The four truth files, indexed for joining, with join integrity checked.

      * truth_dated / truth_tips are keyed by `filename` (one per shipped .nwk)
      * truth_species is keyed by `partition`
      * truth_provenance's canonical rows are keyed by `partition_name` and
        carry the canonical tree's `filename` directly
      * the canonical partitions of truth_provenance are exactly the represented
        partitions of truth_species
    """

    def __init__(self, measured=None):
        for name, path in TRUTH_FILES.items():
            if not os.path.exists(path):
                die("truth file not found: %s\n"
                    "                All four of %s must exist before a rebuild."
                    % (rel(path), ", ".join(os.path.basename(p)
                                            for p in TRUTH_FILES.values())))

        # The measured columns come from refresh_measured_truth() when it has
        # run (always, in main); read from disk otherwise.
        measured = measured or {}
        self.dated_rows = measured.get("dated") or read_csv(TRUTH_FILES["dated"])
        self.tips_rows = measured.get("tips") or read_csv(TRUTH_FILES["tips"])
        self.species_rows = read_csv(TRUTH_FILES["species"])
        self.prov_rows = read_csv(TRUTH_FILES["provenance"])

        for label, rows_ in (("truth_dated", self.dated_rows),
                             ("truth_tips", self.tips_rows),
                             ("truth_species", self.species_rows),
                             ("truth_provenance", self.prov_rows)):
            if not rows_:
                die("%s.csv has no data rows -- a rebuild from an empty truth "
                    "file would silently blank every ledger." % label)

        self.dated = {r["filename"]: r for r in self.dated_rows}
        self.tips = {r["filename"]: r for r in self.tips_rows}
        self.species = OrderedDict((r["partition"], r) for r in self.species_rows)

        self.canon_prov = OrderedDict(
            (r["partition_name"], r) for r in self.prov_rows
            if truthy(r.get("is_canonical")))
        self.alt_prov = [r for r in self.prov_rows
                         if not truthy(r.get("is_canonical"))]

        if set(self.dated) != set(self.tips):
            die("truth_dated.csv and truth_tips.csv disagree on the set of "
                "shipped .nwk files (%d vs %d rows, %d in only one)."
                % (len(self.dated), len(self.tips),
                   len(set(self.dated) ^ set(self.tips))))

        represented = {p for p, r in self.species.items()
                       if not truthy(r.get("is_unrepresented"))}
        if represented != set(self.canon_prov):
            die("the represented partitions of truth_species.csv and the "
                "canonical rows of truth_provenance.csv are not the same set.\n"
                "                represented-only: %s\n"
                "                canonical-only:   %s"
                % (sorted(represented - set(self.canon_prov)),
                   sorted(set(self.canon_prov) - represented)))

        # The canonical tree is named by truth_provenance directly; the
        # group-stem convention is cross-checked, never assumed.
        self.canon_file = OrderedDict()
        for part, pv in self.canon_prov.items():
            fname = (pv.get("filename") or "").strip()
            if fname not in self.tips:
                die("partition %r names canonical tree %r, which is not among "
                    "the shipped trees of truth_tips.csv." % (part, fname))
            self.canon_file[part] = fname

        self.eukaryotic = [p for p, r in self.species.items()
                           if not truthy(r.get("is_prokaryotic"))]
        euk = set(self.eukaryotic)
        self.canon_euk = [p for p in self.canon_prov if p in euk]
        self.canon_prok = [p for p in self.canon_prov if p not in euk]

    # -- per-partition accessors ------------------------------------------
    def group(self, partition):
        return (self.species[partition].get("group") or "").strip()

    def labels_unique(self, partition):
        """The adopted tip convention: distinct standardised labels in the
        shipped Newick file (truth_tips.csv `tips_adopted`)."""
        return as_int(self.tips[self.canon_file[partition]]["tips_adopted"], 0)

    def labels_raw(self, partition):
        return as_int(self.tips[self.canon_file[partition]]["raw_leaves"], 0)

    def is_geometric_chronogram(self, partition):
        """The file's own verdict: root-to-tip CV within tolerance."""
        return truthy(self.dated[self.canon_file[partition]]["is_chronogram"])

    def time_scale(self, partition):
        return time_scales().get(self.canon_file[partition], "")

    def is_chronogram(self, partition):
        """C1d: an ultrametric file whose branch lengths are ABSOLUTE time
        (decision D14).  An ultrametric tree in relative time is reported as
        its own class and does not count as dated."""
        return (self.is_geometric_chronogram(partition)
                and self.time_scale(partition) == "absolute")

    def roottotip_cv(self, partition):
        return as_float(self.dated[self.canon_file[partition]]["roottotip_cv"])

    def archival_class(self, partition):
        return (self.canon_prov[partition].get("archival_class") or "").strip()

    def interval_recoverable(self, partition):
        """MEASURED uncertainty flag (S02, rung C1b): the deposit of the
        canonical tree carries a divergence-time interval for its nodes, of a
        kind that counts (Data/crown_intervals.csv, written by
        scripts/crown_intervals.py from the deposits).  No ledger column
        enters: archival_class and divergence_ci_available are prose."""
        r = crown_intervals().get(partition)
        return bool(r) and r.get("counts_for_c1b") == "yes"

    def interval_kind(self, partition):
        r = crown_intervals().get(partition)
        return r["kind"] if r and r.get("counts_for_c1b") == "yes" else ""

    def posterior_n(self, partition):
        """Number of dated trees behind the interval (samples only)."""
        r = crown_intervals().get(partition)
        if not r or r.get("counts_for_c1b") != "yes" \
                or as_int(r.get("n_trees"), 0) <= 1:
            return ""
        return str(as_int(r["n_trees"]))

    def described(self, partition):
        return as_int(self.species[partition].get("described"), 0)

    def is_unrepresented(self, partition):
        return truthy(self.species[partition].get("is_unrepresented"))

    def is_prokaryotic(self, partition):
        return truthy(self.species[partition].get("is_prokaryotic"))


# ==========================================================================
# THE NAME-RESOLUTION LAYER
# ==========================================================================

class Names(object):
    """The deposit's own name resolution, used for BOTH sides of every
    intersection.

      dictionary.csv    numeric tip id -> standardised name.  The shipped trees
                        carry numeric ids, so this is what turns a shipped tree
                        into a set of species labels.
      full_mapping.csv  (group, original label) -> standardised name.  Raw
                        source trees carry the ORIGINAL labels, and several
                        source deposits decorate them with suffixes that the
                        standardisation strips.  Resolving raw labels through
                        the deposit's own mapping -- rather than through a
                        bespoke suffix heuristic -- is what makes ONE formula
                        work for every partition.
    """

    def __init__(self, warnings):
        for path in (DICTIONARY, FULL_MAPPING):
            if not os.path.exists(path):
                die("missing %s -- the name-resolution layer is required: the "
                    "shipped trees carry numeric ids, and without it no "
                    "numerator can be formed." % rel(path))
        self.by_id = {}
        with open(DICTIONARY, newline="", encoding="utf-8",
                  errors="replace") as fh:
            for row in csv.reader(fh):
                if len(row) < 2 or row[0] == "id":
                    continue
                self.by_id[row[0]] = row[1]
        if not self.by_id:
            die("%s parsed to zero id -> name pairs." % rel(DICTIONARY))

        self.by_group = defaultdict(dict)
        with open(FULL_MAPPING, newline="", encoding="utf-8",
                  errors="replace") as fh:
            for row in csv.DictReader(fh):
                g = row.get("group")
                o = row.get("original_label")
                s = row.get("standardized_name")
                if g and o and s:
                    self.by_group[g][o] = s
        if not self.by_group:
            warnings.append(
                "%s parsed to zero (group, original label) pairs; raw source "
                "trees whose labels differ from the standardised form cannot "
                "be resolved and will be rejected as sources"
                % rel(FULL_MAPPING))

    def resolve_ids(self, labels):
        """Shipped-tree numeric ids -> standardised names."""
        out, unmapped = set(), 0
        for lab in labels:
            name = self.by_id.get(lab)
            if name is None:
                unmapped += 1
                name = lab
            out.add(name)
        return out, unmapped

    def resolve_originals(self, labels, group):
        """Raw-source original labels -> standardised names, via the deposit's
        own mapping for that group.  Labels absent from the mapping pass
        through unchanged, which is correct for the (majority) sources whose
        labels are already in standardised form."""
        m = self.by_group.get(group, {})
        return set(m.get(lab, lab) for lab in labels)


class ShippedTrees(object):
    """Label sets and re-measured ultrametricity for the shipped trees."""

    def __init__(self, truth, names, warnings, verify):
        self.truth = truth
        self.names = names
        self.warnings = warnings
        self.labels = {}          # partition -> set(standardised name)
        self.ids = {}             # partition -> {standardised name -> id}
        self.cv_measured = {}     # partition -> re-measured root-to-tip CV
        # partition -> {counting unit -> display name}: CoL accepted names for a
        # eukaryotic partition (S05), GTDB species clusters (labels) for a
        # prokaryotic one (D12)
        self.c0_names = {}
        self.c0_unresolved = {}
        # partition -> {standardized id -> counting unit} for the tips that
        # resolve; C1 is the set of units of the molecular tips (A01)
        self.unit_of = {}
        # partition -> set of C1 counting units, filled by build_basis
        self.c1_names = {}
        self.col_rows = defaultdict(dict)
        if os.path.exists(COL_RESOLUTION):
            for r in read_csv(COL_RESOLUTION):
                self.col_rows[r["partition"]][r["standardized_id"]] = r
        if not os.path.isdir(STD_TREES):
            die("missing shipped-tree directory %s." % rel(STD_TREES))
        for part in truth.canon_prov:
            self._load(part, verify)

    def _load(self, part, verify):
        fname = self.truth.canon_file[part]
        path = os.path.join(STD_TREES, fname)
        if not os.path.exists(path):
            die("canonical tree %s is named by truth_provenance.csv but is not "
                "deposited at %s." % (fname, rel(path)))
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        raw_labels, depths = tree_facts(text)

        names, unmapped = self.names.resolve_ids(raw_labels)
        self.labels[part] = names
        # id lookup for molecular_tips.csv, deterministic on ties
        first = {}
        for lab in raw_labels:
            nm = self.names.by_id.get(lab, lab)
            if nm not in first or lab < first[nm]:
                first[nm] = lab
        self.ids[part] = first

        if self.truth.is_prokaryotic(part):
            self.c0_names[part] = {nm: nm for nm in names}
            self.c0_unresolved[part] = 0
            self.unit_of[part] = {sid: nm for nm, sid in first.items()}
        else:
            rows = self.col_rows.get(part, {})
            tip_ids = set(first.values())
            if set(rows) != tip_ids:
                die("%s: %s has %d rows for %s but the shipped tree has %d "
                    "tips (%d without a row); run scripts/match_col.py"
                    % (fname, rel(COL_RESOLUTION), len(rows), part,
                       len(tip_ids), len(tip_ids - set(rows))))
            acc = {}
            unit_of = {}
            unres = 0
            for sid in tip_ids:
                r = rows[sid]
                if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
                    acc[r["col_id"]] = r["col_name"]
                    unit_of[sid] = r["col_id"]
                else:
                    unres += 1
            self.c0_names[part] = acc
            self.unit_of[part] = unit_of
            self.c0_unresolved[part] = unres

        if unmapped:
            self.warnings.append(
                "%s: %d of %d shipped tip labels are absent from "
                "dictionary.csv and were carried through unresolved"
                % (part, unmapped, len(raw_labels)))

        # Cross-check the adopted tip convention against the file itself.
        adopted = self.truth.labels_unique(part)
        if len(names) != adopted:
            self.warnings.append(
                "%s: the shipped tree resolves to %d distinct standardised "
                "labels but truth_tips.csv adopts %d. Coverage is computed on "
                "the TRUTH FILE's value; investigate before publishing."
                % (part, len(names), adopted))

        if verify:
            cv = coefficient_of_variation(depths)
            self.cv_measured[part] = cv
            recorded = self.truth.roottotip_cv(part)
            if cv is not None and recorded is not None \
                    and abs(cv - recorded) > CV_RECHECK_TOL:
                self.warnings.append(
                    "%s: re-measured root-to-tip CV %.6f differs from "
                    "truth_dated.csv's %.6f by more than the %.0e recheck "
                    "tolerance. GROUND TRUTH MAY HAVE ROTTED -- do not publish "
                    "until this is resolved."
                    % (part, cv, recorded, CV_RECHECK_TOL))


# ==========================================================================
# AXIS A -- MOLECULAR SUPPORT OF PLACEMENT, MEASURED
# ==========================================================================
#
# f_mol comes from ONE formula.  What varies between partitions is only which
# EVIDENCE the deposit carries, and that is read out of truth_provenance.csv's
# own prose, never out of a list of partition names.  The resolution order is:
#
#   1 SOURCE_TREE_INTERSECTION.  The ground truth's molecular_tip_fraction
#     sentence names a deposited raw tree.  We parse it, resolve its labels
#     through the deposit's own mapping, intersect with the shipped tree, and
#     use the MEASURED intersection.  If the sentence also records a numerator
#     N, the measurement is checked against it and the MEASUREMENT WINS -- a
#     disagreement is warned about, never silently absorbed.  (This is the
#     documented Seed-plant correction: the recorded ratio was taken on RAW
#     LEAVES while the atlas counts DISTINCT LABELS.)
#     If the named path is not a usable tree, the directories it points at are
#     scanned, size-capped, and the candidate whose measured intersection best
#     reproduces N is accepted.
#
#   2 PUBLISHED_FRACTION.  No deposited source tree, but the sentence records a
#     fraction N/M from the source publication.  f_mol = N / labels_unique --
#     i.e. rebased onto the basis the atlas actually counts, not onto the
#     publication's own denominator.  Flagged: an upper bound.
#
#   3 ASSERTED_NO_IMPUTATION.  imputation_status == 'none': the shipped tree is
#     the molecular analysis, so it is its own source and f_mol = 1 by the same
#     formula, not by exception.
#
#   4 UNRESOLVED.  None of the above; the ledger scalar is carried through and
#     the partition is named in the disclosures.
#
# CONTESTED SOURCES.  When the ground truth asserts f_mol = 1 with
# imputation_status = 'none' but a deposited tree it names implies less, the
# point estimate KEEPS the asserted 1 (the ground truth's ruling) and the
# implied value is published as the lower end of the band, f_low, plus its own
# sensitivity line.  Silently picking either number would be a decision this
# script has no standing to make.

RAW_PATH_RE = re.compile(r"(?:raw|molecular_sources)/[A-Za-z0-9_./()+\-]+")
RATIO_RE = re.compile(r"(?<![\d.])(\d[\d,]{1,9})\s*/\s*(\d[\d,]{1,9})(?![\d.])")
FRACTION_RE = re.compile(r"(?<![\d.])(0?\.\d+|1\.0+)(?![\d])")
SENTENCE_SPLIT_RE = re.compile(r"(?<=[.;])\s+")


class MolecularEvidence(object):
    """Resolves f_mol, its basis, its band and its tip set for one partition."""

    __slots__ = ("partition", "f_mol", "f_basis", "f_low", "f_high",
                 "source_file", "tip_names", "recorded_numerator",
                 "contested_source", "contested_f", "notes",
                 "measured_intersection", "exact_label_intersection")

    def __init__(self, partition):
        self.partition = partition
        self.f_mol = None
        self.f_basis = "unresolved"
        self.f_low = None
        self.f_high = None
        self.source_file = ""
        self.tip_names = None          # set, only when measured from a tree
        self.recorded_numerator = None
        self.contested_source = ""
        self.contested_f = None
        self.notes = []
        self.measured_intersection = None       # under the adopted rule
        self.exact_label_intersection = None    # under naive exact matching


def molecular_sentences(note):
    """The sentences of a curatorial note that are ABOUT molecular_tip_fraction.

    Scoping the scrape to these sentences is what stops the parser picking up
    the divergence-CI evidence files, the rejected alternatives and the
    citation-verification asides that also live in these notes.
    """
    return [s for s in SENTENCE_SPLIT_RE.split(note or "")
            if "molecular_tip_fraction" in s]


def pruned_ids(filename):
    out = set()
    for r in read_csv(PRUNED_TIPS, required=False):
        if (r.get("filename") or "").strip() == filename:
            out.add((r.get("standardized_id") or "").strip())
    return out


def resolve_molecular(truth, shipped, names, warnings):
    """Per-partition molecular evidence for every canonical partition."""
    out = OrderedDict()
    evidence = defaultdict(list)
    if os.path.exists(MOL_TIP_EVIDENCE):
        for r in read_csv(MOL_TIP_EVIDENCE):
            evidence[(r.get("partition") or "").strip()].append(r)
    for part in truth.canon_prov:
        ev = MolecularEvidence(part)
        pv = truth.canon_prov[part]
        note = pv.get("conflict_note") or ""
        sents = molecular_sentences(note)
        blob = " ".join(sents)
        shipped_names = shipped.labels[part]
        n_shipped = truth.labels_unique(part)
        ledger_f = as_float(pv.get("molecular_tip_fraction"))
        imputation = (pv.get("imputation_status") or "").strip()

        # ---- the recorded numerator, used as a CHECK on the measurement ----
        for m in RATIO_RE.finditer(blob):
            n = as_int(m.group(1))
            d = as_int(m.group(2))
            if n is not None and d is not None and 0 < n <= d:
                ev.recorded_numerator = n
                break

        # ---- 1. a deposited source tree, measured ---------------------------
        named = []
        for s in sents:
            for p in RAW_PATH_RE.findall(s):
                if p not in named:
                    named.append(p)
        candidates = [os.path.join(ROOT, "Data", p) for p in named
                      if p.lower().endswith(SOURCE_TREE_EXT)]
        if named and not candidates:
            # The ground truth points at a directory but the named artefact is
            # not a tree (e.g. a species list).  Look in that directory.
            seen = set()
            for p in named:
                d = os.path.dirname(os.path.join(ROOT, "Data", p))
                if d in seen or not os.path.isdir(d):
                    continue
                seen.add(d)
                for fn in sorted(os.listdir(d)):
                    if fn.lower().endswith(SOURCE_TREE_EXT):
                        candidates.append(os.path.join(d, fn))

        excluded = frozenset(names.by_id.get(i, i)
                             for i in pruned_ids(truth.canon_file[part]))
        best = _best_source(candidates, shipped_names, names,
                            truth.group(part), ev.recorded_numerator, ev,
                            excluded)

        if best is not None:
            inter, path = best
            f = inter / float(n_shipped) if n_shipped else None
            ev.measured_intersection = inter
            ev.exact_label_intersection = len(
                _source_sets(path, names, truth.group(part))[1] & shipped_names)
            if ev.exact_label_intersection != inter:
                warnings.append(
                    "%s: the adopted rule (both label sets resolved through the "
                    "deposit's own name mapping) measures %d molecular tips "
                    "against %s; naive EXACT-LABEL matching measures %d. The "
                    "adopted rule is used because exact matching is a "
                    "measurement artefact -- it scores a decorated-label source "
                    "at or near zero. Both figures are published in "
                    "headline_numbers.json -> molecular_evidence."
                    % (part, inter, rel(path), ev.exact_label_intersection))
            asserted_one = (ledger_f is not None and ledger_f >= 1.0
                            and imputation == "none")
            if asserted_one and f is not None and f < 1.0:
                # CONTESTED: the ground truth ruled for 1.00 and flagged the
                # disagreement.  Keep the ruling; publish the band.
                ev.contested_source = rel(path)
                ev.contested_f = f
                ev.notes.append(
                    "a deposited tree named in the ground truth implies "
                    "f_mol=%.6f against the asserted 1.0; the asserted value "
                    "is kept and the implied value published as f_low"
                    % f)
                best = None
            else:
                ev.f_mol = f
                ev.f_basis = ("source_alignment_intersection"
                              if path.endswith(TAXA_LIST_EXT)
                              else "source_tree_intersection")
                ev.source_file = rel(path)
                src_keys = _species_keys(_source_names(path, names,
                                                       truth.group(part)))
                ev.tip_names = {nm for nm in shipped_names
                                if (_species_keys([nm]) & src_keys)}
                if ev.recorded_numerator is not None \
                        and ev.recorded_numerator != inter:
                    warnings.append(
                        "%s: MEASURED intersection %d against %s differs from "
                        "the numerator %d recorded in truth_provenance.csv. "
                        "The MEASUREMENT is used (the recorded ratio is on a "
                        "different tip basis); truth_provenance.csv's "
                        "molecular_tip_fraction for this partition should be "
                        "updated with a note."
                        % (part, inter, rel(path), ev.recorded_numerator))

        # ---- 1b. per-tip evidence -------------------------------------------
        # Data/molecular_tip_evidence.csv lists, tip by tip, the evidence that
        # a shipped tip is in the molecular analysis behind the tree.  A
        # source_tree_label row is checked against the named deposited tree
        # whenever that tree is present; a row that fails is not counted.
        rows_ev = evidence.get(part) or []
        if ev.f_mol is None and rows_ev and n_shipped:
            names_ok, missing_src = set(), set()
            for r in rows_ev:
                sid = (r.get("standardized_id") or "").strip()
                kind = (r.get("evidence_kind") or "").strip()
                nm = names.by_id.get(sid)
                if kind not in MOL_EVIDENCE_KINDS:
                    die("%s: unknown evidence_kind %r in %s"
                        % (part, kind, rel(MOL_TIP_EVIDENCE)))
                if nm is None or nm not in shipped_names:
                    warnings.append("%s: evidence row for id %s does not "
                                    "resolve to a shipped tip; not counted"
                                    % (part, sid))
                    continue
                if kind == "source_tree_label":
                    src = os.path.join(ROOT, "Data",
                                       (r.get("source_file") or "").strip())
                    if os.path.exists(src):
                        if (r.get("source_label") or "").strip() not in \
                                _source_sets(src, names, truth.group(part))[1]:
                            warnings.append(
                                "%s: id %s cites label %r, which is not in %s; "
                                "not counted" % (part, sid, r.get("source_label"),
                                                 rel(src)))
                            continue
                    elif src not in missing_src:
                        missing_src.add(src)
                        warnings.append(
                            "%s: %s is not present, so the source labels cited "
                            "in %s were not re-checked on this run"
                            % (part, rel(src), rel(MOL_TIP_EVIDENCE)))
                names_ok.add(nm)
            ev.f_mol = len(names_ok) / float(n_shipped)
            ev.f_basis = "per_tip_evidence"
            ev.source_file = rel(MOL_TIP_EVIDENCE)
            ev.tip_names = names_ok
            ev.measured_intersection = len(names_ok)

        # ---- 2. a published fraction ----------------------------------------
        if ev.f_mol is None and ev.recorded_numerator is not None and n_shipped:
            ev.f_mol = min(1.0, ev.recorded_numerator / float(n_shipped))
            ev.f_basis = "published_fraction"
            ev.notes.append(
                "no molecular source tree is deposited; the fraction is the "
                "source publication's, rebased from its own denominator onto "
                "the %d distinct labels the atlas counts. AN UPPER BOUND."
                % n_shipped)

        # ---- 3. asserted absence of imputation -------------------------------
        if ev.f_mol is None and imputation == "none":
            ev.f_mol = 1.0
            ev.f_basis = "asserted_no_imputation"
            ev.notes.append(
                "imputation_status='none': the shipped tree IS the molecular "
                "analysis, so it is its own source and scores 1 by the same "
                "formula, not by exception.")

        # ---- 4. unresolved ---------------------------------------------------
        if ev.f_mol is None:
            ev.f_mol = ledger_f
            ev.f_basis = "unresolved"
            if ev.f_mol is None:
                die("partition %r has no molecular_tip_fraction and no "
                    "resolvable molecular evidence; C1 cannot be formed."
                    % part)
            ev.notes.append(
                "neither a deposited source tree nor a recorded ratio nor an "
                "asserted absence of imputation; the pre-audit ledger scalar "
                "is carried through and the partition is disclosed as "
                "unresolved.")

        # ---- the band -------------------------------------------------------
        lo = hi = ev.f_mol
        if ev.contested_f is not None:
            lo = min(lo, ev.contested_f)
            hi = max(hi, ev.contested_f)
        # An explicitly UNRESOLVED ledger disagreement publishes its competing
        # values as the band.  Keyed on the ground truth's own vocabulary, so
        # rejected alternatives ("REPLACES 0.67") are correctly excluded.
        for s in sents:
            if "UNRESOLVED" not in s:
                continue
            vals = [as_float(v) for v in FRACTION_RE.findall(s)]
            vals = [v for v in vals if v is not None and 0.0 <= v <= 1.0]
            for v in vals:
                lo = min(lo, v)
                hi = max(hi, v)
        ev.f_low, ev.f_high = lo, hi
        if hi > lo:
            ev.notes.append("BAND: f_mol in [%.6f, %.6f]" % (lo, hi))

        if ev.f_mol < 0.0 or ev.f_mol > 1.0:
            die("partition %r resolved to f_mol=%r, outside [0,1]."
                % (part, ev.f_mol))
        out[part] = ev
    return out


_SOURCE_CACHE = {}


def _source_sets(path, names, group):
    """(resolved labels, verbatim labels) of a deposited source tree.

    RESOLVED is the adopted rule: raw labels are passed through the deposit's
    OWN (group, original label) -> standardised name mapping before the
    intersection, exactly as the shipped tree's ids are.  Resolving BOTH sides
    with the same layer is what makes one formula work for every partition:
    several source deposits decorate their labels (collector numbers, voucher
    suffixes), and naive exact matching silently scores such a source ZERO --
    a measurement artefact, not a biological fact.

    VERBATIM is kept only as a published diagnostic, so a reader can reproduce
    the naive figure and see exactly what the resolution step changed.
    """
    key = (path, group)
    if key in _SOURCE_CACHE:
        return _SOURCE_CACHE[key]
    empty = (frozenset(), frozenset())
    if path.endswith(TAXA_LIST_EXT):
        with open(path, encoding="utf-8") as fh:
            labels = [ln.strip() for ln in fh
                      if ln.strip() and not ln.startswith("#")]
        _SOURCE_CACHE[key] = (frozenset(names.resolve_originals(labels, group)),
                              frozenset(labels))
        return _SOURCE_CACHE[key]
    text = read_first_tree(path, SOURCE_SCAN_MAX_BYTES)
    if not text:
        _SOURCE_CACHE[key] = empty
        return empty
    try:
        labels, _ = tree_facts(text)
    except Exception:                                     # noqa: BLE001
        _SOURCE_CACHE[key] = empty
        return empty
    _SOURCE_CACHE[key] = (frozenset(names.resolve_originals(labels, group)),
                          frozenset(labels))
    return _SOURCE_CACHE[key]


def _source_names(path, names, group):
    return _source_sets(path, names, group)[0]


def _species_keys(labels):
    """Decision D4: a tip stands for its species.  Each label is reduced to its
    species key (scripts/canonical_trees.py species_key); a label that is not a
    species label keeps itself as its key, so it matches only itself."""
    out = set()
    for nm in labels:
        k, ok = SPECIES_KEY(nm, nm)
        out.add(k if ok else "label:" + nm)
    return out


def _best_source(candidates, shipped_names, names, group, recorded, ev,
                 excluded=frozenset()):
    """Pick the deposited tree that best evidences the molecular subset.

    ACCEPTANCE (both must hold):
      * the source's labels are (near) a SUBSET of the shipped tree's -- a
        source that is mostly outside the shipped tree is a different taxon
        sample, not the analysis behind this tree;
      * where the ground truth recorded a numerator, the measured intersection
        reproduces it within SOURCE_NUMERATOR_TOL.
    Ties break on the smaller |measured - recorded|, then on the file path, so
    the choice is deterministic and a rerun cannot flip it.
    """
    scored = []
    for path in candidates:
        if not os.path.exists(path):
            continue
        # labels pruned from the canonical tree (D7) are known to lie outside
        # it; they are not evidence that the source is a different sample
        src = _species_keys(_source_names(path, names, group) - excluded)
        if not src:
            continue
        shipped_keys = _species_keys(shipped_names)
        inter = len(src & shipped_keys)
        if inter == 0:
            ev.notes.append(
                "%s shares no resolved label with the shipped tree and was "
                "rejected as a source" % rel(path))
            continue
        outside = len(src - shipped_keys) / float(len(src))
        if outside > SOURCE_NONSUBSET_TOL:
            ev.notes.append(
                "%s has %.1f%% of its labels outside the shipped tree "
                "(tolerance %.0f%%) and was rejected as a source"
                % (rel(path), 100.0 * outside, 100.0 * SOURCE_NONSUBSET_TOL))
            continue
        if recorded is not None:
            dev = abs(inter - recorded) / float(recorded)
            if dev > SOURCE_NUMERATOR_TOL:
                ev.notes.append(
                    "%s implies %d molecular tips against the %d the ground "
                    "truth records (%.1f%% off, tolerance %.0f%%); rejected"
                    % (rel(path), inter, recorded, 100.0 * dev,
                       100.0 * SOURCE_NUMERATOR_TOL))
                continue
            scored.append((dev, path, inter))
        else:
            scored.append((0.0, path, inter))
    if not scored:
        return None
    scored.sort(key=lambda t: (t[0], t[1]))
    dev, path, inter = scored[0]
    return inter, path


# ==========================================================================
# THE PRE-AUDIT LEDGER -- read for the NOMINAL annotation only
# ==========================================================================

class PreAuditLedger(object):
    """Data/tree_provenance.csv, demoted.

    ONLY `inference`, `otl_pipeline_involvement` and `tree_class` are read, and
    `tree_class` only to reprint S5's descriptive column.  Every quantity that
    this ledger also carries -- dated, ntips, molecular_tip_fraction -- is
    ignored in favour of the measured truth files.  A paper whose thesis is
    "ledger flags must be measured" cannot certify a reporting axis on
    uncertified ledger flags, which is exactly why topology_assembly is a
    NOMINAL ANNOTATION and never a rung.
    """

    def __init__(self, truth, warnings):
        rows = read_csv(TREE_PROVENANCE, required=False)
        self.by_partition = {}
        for r in rows:
            if truthy(r.get("is_partition_canonical")):
                self.by_partition[r.get("partition_name")] = r
        if not self.by_partition:
            warnings.append(
                "%s carries no partition-canonical rows; topology_assembly "
                "falls back to '%s' for every partition and the Figure 2 "
                "marker shapes lose their meaning"
                % (rel(TREE_PROVENANCE), ASSEMBLY_PRIMARY))
        missing = [p for p in truth.canon_prov if p not in self.by_partition]
        if missing:
            warnings.append(
                "no pre-audit ledger row for %d canonical partition(s): %s; "
                "topology_assembly defaults to '%s' there"
                % (len(missing), ", ".join(sorted(missing)), ASSEMBLY_PRIMARY))

    def triple(self, partition):
        r = self.by_partition.get(partition) or {}
        return ((r.get("otl_pipeline_involvement") or "").strip(),
                (r.get("inference") or "").strip(),
                (r.get("tree_class") or "").strip())

    def assembly(self, partition):
        otl, inference, _ = self.triple(partition)
        return classify_assembly(otl, inference)

    def placement_evidence(self, partition):
        """The raw triple, shipped verbatim so any reader can overturn the
        classifier in four lines."""
        otl, inference, tree_class = self.triple(partition)
        return "otl=%s|inference=%s|tree_class=%s" % (otl or "-",
                                                      inference or "-",
                                                      tree_class or "-")

    def tree_class(self, partition):
        return self.triple(partition)[2].replace("_", " ").strip()

    def notes(self, partition):
        r = self.by_partition.get(partition) or {}
        return (r.get("notes") or "").strip()


# ==========================================================================
# PASSTHROUGH resolution (fields the truth files deliberately do not carry)
# ==========================================================================

class Passthrough(object):
    """estimated_total and confidence, by a documented precedence chain.

    Neither is derivable from the four truth files, so both are reported as
    PASSTHROUGH.  estimated_total is a MAJORITY VOTE across the three ledgers
    that carry it -- deliberately a vote, not a precedence chain, because that
    mirrors the rule truth_species.csv itself used to resolve `described`.
    """

    def __init__(self, truth, warnings):
        self.truth = truth
        self.warnings = warnings
        self._yaml = self._load_yaml_module()
        self._info = {}
        self._total_cache = {}
        self._slug = {}
        self._resolve_info_dirs()
        self.data_estimates = {r.get("group"): r for r in
                               read_csv(DATA_ESTIMATES, required=False)}
        self.old_summary = {r.get("partition"): r for r in
                            read_csv(OLD_PARTITION_SUMMARY, required=False)}

    @staticmethod
    def _load_yaml_module():
        try:
            import yaml           # noqa: PLC0415 - optional dependency
            return yaml
        except ImportError:
            return None

    def _resolve_info_dirs(self):
        if not os.path.isdir(PARTITION_INFO_DIR):
            self.warnings.append(
                "curated layer %s not found; estimated_total and confidence "
                "fall back to data_estimates.csv" % rel(PARTITION_INFO_DIR))
            return
        present = {d for d in os.listdir(PARTITION_INFO_DIR)
                   if os.path.exists(os.path.join(PARTITION_INFO_DIR, d,
                                                  "info.yaml"))}
        for part, row in self.truth.species.items():
            for cand in (slugify(part), (row.get("group") or "")):
                if cand in present:
                    self._slug[part] = cand
                    break

    def _yaml_estimate(self, partition):
        if self._yaml is None or partition not in self._slug:
            return {}
        if partition in self._info:
            return self._info[partition]
        path = os.path.join(PARTITION_INFO_DIR, self._slug[partition],
                            "info.yaml")
        try:
            with open(path, encoding="utf-8") as fh:
                doc = self._yaml.safe_load(fh) or {}
        except Exception as exc:                          # noqa: BLE001
            self.warnings.append("could not parse %s (%s)" % (rel(path), exc))
            doc = {}
        self._info[partition] = doc.get("estimate") or {}
        return self._info[partition]

    def estimated_total(self, partition):
        if partition in self._total_cache:
            return self._total_cache[partition]
        est = self._yaml_estimate(partition)
        curated = as_int(est.get("total"))
        candidates = [
            ("curated info.yaml", curated),
            ("data_estimates.csv",
             as_int((self.data_estimates.get(partition) or {})
                    .get("estimated_total"))),
            ("partition_summary.csv",
             as_int((self.old_summary.get(partition) or {})
                    .get("estimated_total"))),
        ]
        votes = defaultdict(list)
        for src, val in candidates:
            if val is not None:
                votes[val].append(src)
        if not votes:
            self.warnings.append(
                "estimated_total unresolved for partition %r; emitted blank"
                % partition)
            self._total_cache[partition] = None
            return None
        best = max(votes, key=lambda v: (len(votes[v]),
                                         1 if v == curated else 0))
        if len(votes) > 1:
            self.warnings.append(
                "estimated_total for %r: took %s (%s) over %s"
                % (partition, commas(best), ",".join(votes[best]),
                   "; ".join("%s=%s" % (",".join(votes[v]), commas(v))
                             for v in sorted(votes) if v != best)))
        self._total_cache[partition] = best
        return best

    def confidence(self, partition):
        v = (self._yaml_estimate(partition) or {}).get("confidence")
        if v:
            return str(v)
        return ((self.data_estimates.get(partition) or {})
                .get("confidence") or "").strip()


# ==========================================================================
# THE BASIS TABLE -- 62 rows, one per partition.  Everything else is a sum
# over columns of THIS table.
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

def c1_count(part, ev, shipped):
    """C1 for one partition, in the same unit as C0 (audit A01).

    C1 is the number of distinct counting units (CoL accepted names; GTDB
    species clusters for prokaryotes) among the tips that are in the molecular
    analysis.  Where the molecular tips are enumerated it is that count
    ('enumerated').  Where f_mol = 1 every tip is molecular and C1 = C0
    ('all_tips').  Otherwise it is the expectation c0_units * f_mol
    ('expectation'), which assumes molecular and non-molecular labels resolve
    to accepted names at the same rate -- an assumption the enumerated
    partitions show to be false, so it is used only where nothing better
    exists.  Returns (units, basis, set of units or None)."""
    c0 = shipped.c0_names.get(part, {})
    if ev.tip_names is not None and ev.source_file:
        idmap = shipped.ids.get(part, {})
        unit_of = shipped.unit_of.get(part, {})
        units = set()
        for nm in ev.tip_names:
            sid = idmap.get(nm)
            if sid is not None and sid in unit_of:
                units.add(unit_of[sid])
        return len(units), "enumerated", units
    if ev.f_mol is not None and ev.f_mol == 1.0:
        return len(c0), "all_tips", set(c0)
    f = ev.f_mol if ev.f_mol is not None else 0.0
    return round(len(c0) * f, NUMERATOR_DP), "expectation", None


def build_basis(truth, shipped, mol, ledger, warnings):
    """The 62-row basis table, as an ordered list of dicts (not yet CSV text).

    THE 36 PARTITIONS WITH NO TREE ARE ROWS, with labels_unique = 0.  That is
    the structural point of the file: a numerator and its denominator are then
    column sums of the SAME table, so a denominator that exists in no file
    cannot be written down.
    """
    rows = []
    for part, sp in truth.species.items():
        unrep = truth.is_unrepresented(part)
        row = OrderedDict([
            ("partition", part),
            ("group", (sp.get("group") or "").strip()),
            ("category", (sp.get("category") or "").strip()),
            ("is_prokaryotic", "yes" if truth.is_prokaryotic(part) else "no"),
            ("described", truth.described(part)),
            ("described_source", (sp.get("described_source") or "").strip()),
            ("is_unrepresented", "yes" if unrep else "no"),
            ("canonical_tree", ""),
            ("labels_raw", 0),
            ("labels_unique", 0),
            ("c0_units", 0),
            ("c0_unresolved_labels", 0),
            ("f_mol", ""),
            ("f_basis", ""),
            ("f_low", ""),
            ("f_high", ""),
            ("c1_units", ""),
            ("c1_basis", ""),
            ("molecular_source_file", ""),
            ("is_chronogram", ""),
            ("time_scale", ""),
            ("roottotip_cv", ""),
            ("interval_recoverable", ""),
            ("interval_kind", ""),
            ("posterior_n", ""),
            ("topology_assembly", ""),
            ("placement_evidence", ""),
        ])
        if not unrep:
            ev = mol[part]
            row["canonical_tree"] = truth.canon_file[part]
            row["labels_raw"] = truth.labels_raw(part)
            row["labels_unique"] = truth.labels_unique(part)
            row["c0_units"] = len(shipped.c0_names[part])
            row["c0_unresolved_labels"] = shipped.c0_unresolved[part]
            row["f_mol"] = fnum(ev.f_mol, FRAC_DP)
            row["f_basis"] = ev.f_basis
            row["f_low"] = fnum(ev.f_low, FRAC_DP)
            row["f_high"] = fnum(ev.f_high, FRAC_DP)
            c1, c1_basis, c1_set = c1_count(part, ev, shipped)
            shipped.c1_names[part] = c1_set
            row["c1_units"] = (str(c1) if isinstance(c1, int)
                               else fnum(c1, NUMERATOR_DP))
            row["c1_basis"] = c1_basis
            row["molecular_source_file"] = ev.source_file
            row["is_chronogram"] = "yes" if truth.is_chronogram(part) else "no"
            row["time_scale"] = (truth.time_scale(part) or "") \
                if truth.is_geometric_chronogram(part) else "not_ultrametric"
            if row["time_scale"] not in ("absolute", "relative",
                                         "not_ultrametric"):
                die("%s is ultrametric but has no absolute/relative row in %s "
                    "(decision D14)" % (truth.canon_file[part],
                                        rel(TIME_SCALE_CHECKS)))
            row["roottotip_cv"] = fnum(truth.roottotip_cv(part), FRAC_DP)
            row["interval_recoverable"] = ("yes" if truth.interval_recoverable(part)
                                           else "no")
            row["interval_kind"] = truth.interval_kind(part)
            row["posterior_n"] = truth.posterior_n(part)
            row["topology_assembly"] = ledger.assembly(part)
            row["placement_evidence"] = ledger.placement_evidence(part)
        rows.append(row)
    rows.sort(key=lambda r: r["partition"])
    return rows


# --------------------------------------------------------------------------
# atlas_coverage() -- the reader shipped in Data Availability.  Every published
# figure is one call of it.  Implemented here so that this script and the
# deposit's reader cannot drift.
# --------------------------------------------------------------------------

def atlas_coverage(basis, predicate=None, weight_col=None,
                   denominator_predicate=None, units_col="c0_units"):
    """Coverage from the basis table.

    basis                  list of dict rows of coverage_basis.csv
    predicate              row -> bool, which partitions contribute a NUMERATOR
    weight_col             column of a weight in [0,1]; None means weight 1
    units_col              the count summed: c0_units (C0) or c1_units (C1)
    denominator_predicate  row -> bool, which partitions contribute to D

    Returns dict(numerator, denominator, percent, n_partitions).
    """
    if predicate is None:
        def predicate(_r):
            return True
    if denominator_predicate is None:
        denominator_predicate = predicate
    num = 0.0
    n = 0
    for r in basis:
        if not predicate(r):
            continue
        w = 1.0 if weight_col is None else (as_float(r.get(weight_col)) or 0.0)
        num += (as_float(r.get(units_col)) or 0.0) * w
        n += 1
    den = sum((as_int(r.get("described"), 0) or 0)
              for r in basis if denominator_predicate(r))
    return {"numerator": num, "denominator": den,
            "percent": pct(num, den), "n_partitions": n}


ATLAS_COVERAGE_PY = r'''"""atlas_coverage -- the whole coverage layer of the phylogenetic atlas.

Every published coverage figure in the paper is ONE call of this function
against Data/coverage_basis.csv.  Numerator and denominator are column sums of
the SAME 62-row file, which is what makes a phantom denominator impossible.

    from atlas_coverage import load_basis, atlas_coverage
    b = load_basis("Data/coverage_basis.csv")
    euk = lambda r: r["is_prokaryotic"] == "no"
    atlas_coverage(b, euk)                                   # C0
    atlas_coverage(b, euk, units_col="c1_units")             # C1
    atlas_coverage(b, lambda r: euk(r) and r["is_chronogram"] == "yes",
                   units_col="c1_units", denominator_predicate=euk)   # C1d
    atlas_coverage(b, lambda r: euk(r) and r["interval_recoverable"] == "yes",
                   units_col="c1_units", denominator_predicate=euk)   # C1b

c1_units is the number of distinct accepted names among the molecular tips
(c1_basis = enumerated), c0_units where every tip is molecular (all_tips), and
the expectation c0_units * f_mol only where neither holds (expectation).
"""
import csv


def load_basis(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(x, default=0.0):
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return default


def atlas_coverage(basis, predicate=None, weight_col=None,
                   denominator_predicate=None, units_col="c0_units"):
    predicate = predicate or (lambda r: True)
    denominator_predicate = denominator_predicate or predicate
    num, n = 0.0, 0
    for r in basis:
        if not predicate(r):
            continue
        w = 1.0 if weight_col is None else _num(r.get(weight_col), 0.0)
        num += _num(r.get(units_col)) * w
        n += 1
    den = sum(_num(r.get("described")) for r in basis
              if denominator_predicate(r))
    return {"numerator": num, "denominator": den, "n_partitions": n,
            "percent": (100.0 * num / den) if den else None}
'''

ATLAS_COVERAGE_R = r'''# atlas_coverage -- the whole coverage layer of the phylogenetic atlas.
#
# Every published coverage figure in the paper is ONE call of this function
# against Data/coverage_basis.csv.  Numerator and denominator are column sums
# of the SAME 62-row file, which is what makes a phantom denominator
# impossible.
#
#   basis <- load_basis("Data/coverage_basis.csv")
#   euk   <- function(d) d$is_prokaryotic == "no"
#   atlas_coverage(basis, euk)                                    # C0
#   atlas_coverage(basis, euk, units_col = "c1_units")            # C1
#   atlas_coverage(basis, function(d) euk(d) & d$is_chronogram == "yes",
#                  units_col = "c1_units", denominator_predicate = euk)  # C1d
#
# c1_units is the number of distinct accepted names among the molecular tips
# (c1_basis = enumerated), c0_units where every tip is molecular (all_tips),
# and the expectation c0_units * f_mol only where neither holds.

load_basis <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

atlas_coverage <- function(basis, predicate = NULL, weight_col = NULL,
                           denominator_predicate = NULL,
                           units_col = "c0_units") {
  if (is.null(predicate)) predicate <- function(d) rep(TRUE, nrow(d))
  if (is.null(denominator_predicate)) denominator_predicate <- predicate
  keep <- predicate(basis)
  keep[is.na(keep)] <- FALSE
  w <- if (is.null(weight_col)) rep(1, nrow(basis)) else
    suppressWarnings(as.numeric(basis[[weight_col]]))
  w[is.na(w)] <- 0
  tips <- suppressWarnings(as.numeric(basis[[units_col]]))
  tips[is.na(tips)] <- 0
  num <- sum(tips[keep] * w[keep])
  dkeep <- denominator_predicate(basis)
  dkeep[is.na(dkeep)] <- FALSE
  desc <- suppressWarnings(as.numeric(basis$described))
  desc[is.na(desc)] <- 0
  den <- sum(desc[dkeep])
  list(numerator = num, denominator = den, n_partitions = sum(keep),
       percent = if (den > 0) 100 * num / den else NA_real_)
}
'''


# ==========================================================================
# THE RUNGS
# ==========================================================================

def rung_weight(basis_row, rung):
    """w(rung) for one basis row.  The rungs are NESTED CONJUNCTS: each is the
    previous one AND one more measured condition, so the weights are
    non-increasing by construction."""
    if as_int(basis_row.get("c0_units"), 0) <= 0:
        return 0.0
    w = 1.0
    if rung in ("C1", "C1d", "C1b"):
        w *= c1_fraction(basis_row)
    if rung in ("C1d", "C1b"):
        w *= 1.0 if basis_row.get("is_chronogram") == "yes" else 0.0
    if rung == "C1b":
        w *= 1.0 if basis_row.get("interval_recoverable") == "yes" else 0.0
    return w


def c1_fraction(basis_row, product=False):
    """c1_units / c0_units, so that c0_units * w is the C1 count exactly.
    product=True gives the superseded definition, the label fraction f_mol
    applied to name units, kept as a sensitivity (audit A01)."""
    c0 = as_int(basis_row.get("c0_units"), 0)
    c1 = as_float(basis_row.get("c1_units"))
    if product or c1 is None or not c0:
        return as_float(basis_row.get("f_mol")) or 0.0
    return c1 / float(c0)


def is_euk(row):
    return row.get("is_prokaryotic") == "no"


def compute_rungs(basis, shipped, truth, denominator=None, keep=None):
    """Every rung, over the eukaryotic partitions of `basis`.

    `keep` optionally restricts which partitions may contribute a NUMERATOR
    (used by the leave-assembly-class-out sensitivity); the denominator is
    unaffected, because dropping a partition from the numerator does not make
    its species stop existing.
    """
    euk = [r for r in basis if is_euk(r)]
    den = denominator if denominator is not None \
        else sum(as_int(r["described"], 0) for r in euk)
    contributing = [r for r in euk
                    if as_int(r.get("c0_units"), 0) > 0
                    and (keep is None or keep(r))]

    out = OrderedDict()
    for rung in RUNGS:
        # Rounded ONCE, here, and every downstream quantity is formed from the
        # rounded value.  That is what makes headline_numbers.json exactly
        # self-consistent: every printed percentage is reproducible from the
        # numerator and denominator printed beside it, to the last digit.
        total = round(sum(as_int(r["c0_units"], 0) * rung_weight(r, rung)
                          for r in contributing), NUMERATOR_DP)
        entry = OrderedDict([
            ("definition", RUNGS[rung]),
            ("numerator_sum", total),
            ("numerator_sum_rounded", int(round(total))),
            ("denominator", den),
            ("percent", pct(total, den)),
            ("n_partitions_contributing",
             sum(1 for r in contributing if rung_weight(r, rung) > 0)),
            ("unit", "species" if rung == "C0" else "tip-equivalents"),
        ])
        if rung != "C0":
            # C1 and below are SETS wherever the molecular tips are enumerated
            # or every tip is molecular (c1_basis enumerated / all_tips), so a
            # union is defined when every contributing partition is one.
            live = [r for r in contributing if rung_weight(r, rung) > 0]
            entry["c1_bases"] = OrderedDict(
                (b, sum(1 for r in live if r.get("c1_basis") == b))
                for b in ("enumerated", "all_tips", "expectation"))
            sets = [shipped.c1_names.get(r["partition"]) for r in live]
            # the union is defined only where each row's c1_units IS its set
            # (not where a sensitivity has rescaled it)
            if live and all(x is not None
                            and as_float(r.get("c1_units")) == float(len(x))
                            for r, x in zip(live, sets)):
                union = set().union(*sets)
                entry["numerator_union"] = len(union)
                entry["percent_union"] = pct(len(union), den)
                entry["shared_labels"] = int(round(total)) - len(union)
        if rung == "C0":
            # C0 is the only rung that is a genuine SET, so it is the only one
            # for which a union is defined.  Both conventions are printed so
            # they can never drift apart in the prose.
            parts = [r["partition"] for r in contributing]
            sets = [set(shipped.c0_names[p]) for p in parts
                    if p in shipped.c0_names]
            union = set().union(*sets) if sets else set()
            entry["numerator_union"] = len(union)
            entry["percent_union"] = pct(len(union), den)
            entry["shared_labels"] = int(round(total)) - len(union)
            entry["headline_convention"] = "union"
        out[rung] = entry
    return out


def shared_label_detail(basis, shipped, limit=None):
    """The labels that occur in more than one eukaryotic canonical partition.

    These are outgroups and homonyms.  They are the difference between the
    per-partition SUM and the UNION, and every derived quantity in the paper
    (dark matter above all) must be formed on ONE of the two conventions
    consistently.
    """
    where = defaultdict(list)
    display = {}
    for r in basis:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        p = r["partition"]
        for lab, disp in shipped.c0_names.get(p, {}).items():
            where[lab].append(p)
            display[lab] = disp
    dupes = OrderedDict()
    for lab in sorted(where):
        if len(where[lab]) > 1:
            dupes[display.get(lab, lab)] = sorted(where[lab])
    items = list(dupes.items())
    if limit is not None:
        items = items[:limit]
    return dupes, OrderedDict(items)


# ==========================================================================
# SENSITIVITIES -- all four accompany the headlines and are MANDATORY
# ==========================================================================

def band_sensitivity(basis, shipped, truth, warnings):
    """§8.1 / §8.4 -- every partition whose f_mol is banded, one line each.

    A partition is banded either because the ground truth records an explicitly
    UNRESOLVED disagreement (the deposit's largest open uncertainty) or because
    a deposited tree contests an asserted 1.0.  Both are computed the same way:
    substitute the band end and recompute every rung.
    """
    out = OrderedDict()
    for row in basis:
        if not is_euk(row) or as_int(row.get("c0_units"), 0) <= 0:
            continue
        f = as_float(row.get("f_mol"))
        lo = as_float(row.get("f_low"))
        hi = as_float(row.get("f_high"))
        if f is None or lo is None or hi is None or hi <= lo:
            continue
        entry = OrderedDict([
            ("partition", row["partition"]),
            ("f_point", f),
            ("f_low", lo),
            ("f_high", hi),
            ("f_basis", row.get("f_basis")),
            ("is_chronogram", row.get("is_chronogram")),
        ])
        for tag, val in (("at_f_low", lo), ("at_f_high", hi)):
            alt = [dict(r) for r in basis]
            for r in alt:
                if r["partition"] == row["partition"]:
                    r["f_mol"] = fnum(val, FRAC_DP)
                    # C1 is a count; the band end moves it in proportion
                    c1 = as_float(r.get("c1_units"))
                    if c1 is not None and f:
                        r["c1_units"] = fnum(c1 * val / f, NUMERATOR_DP)
            entry[tag] = _rung_percents(compute_rungs(alt, shipped, truth))
        # Does the band propagate past C1?  If the partition is not a
        # chronogram it CANNOT, and that is the single most reassuring fact
        # about the deposit's largest open uncertainty -- so it is stated.
        entry["propagates_to_C1d"] = (row.get("is_chronogram") == "yes")
        entry["note"] = (
            "the band moves C1; it does NOT move C1d or C1b, because this "
            "partition's canonical tree is not a measured chronogram"
            if row.get("is_chronogram") != "yes"
            else "the band moves C1, C1d and (if the interval is recoverable) C1b")
        out[row["partition"]] = entry
    if not out:
        warnings.append(
            "no banded partition was found; the f_mol band sensitivity is "
            "empty. If the ground truth still records an UNRESOLVED molecular "
            "fraction, the band parser has stopped matching it.")
    return out


def _rung_percents(rungs):
    """A sensitivity snapshot of the ladder.

    REACHABILITY.  Each rung publishes the EXACT (fractional) numerator, the
    rounded one and the denominator, so `percent` is reproducible from the pair
    printed beside it.  Printing only the rounded integer would leave every
    sensitivity percentage unreachable from anything the file carries -- the
    same defect as the retired 23.7 % headline, one level down.
    """
    o = OrderedDict()
    for k, v in rungs.items():
        o[k] = OrderedDict([("numerator", v["numerator_sum_rounded"]),
                            ("numerator_sum", v["numerator_sum"]),
                            ("denominator", v["denominator"]),
                            ("percent", v["percent"])])
        if "percent_union" in v:
            o[k]["numerator_union"] = v["numerator_union"]
            o[k]["percent_union"] = v["percent_union"]
    return o


GBIF_UNMATCHED = os.path.join(ROOT, "Data", "external", "gbif_unmatched.csv")
GBIF_RESOLVED_MATCH = ("EXACT", "FUZZY")
GBIF_RESOLVED_STATUS = ("ACCEPTED", "SYNONYM", "HETEROTYPIC_SYNONYM",
                        "HOMOTYPIC_SYNONYM", "PROPARTE_SYNONYM", "DOUBTFUL")


def col_canonical(col_name):
    """'Genus epithet' of a CoL species-rank name with authorship."""
    return " ".join((col_name or "").split()[:2])


def c0_gbif_sensitivity(basis, shipped):
    """Audit A10.  C0 counts only names in the COL26.9 base release, while each
    denominator comes from its own authority.  Here the tips that release
    leaves 'unmatched' count too when the GBIF backbone resolves their query
    to a species (Data/external/gbif_unmatched.csv, scripts/gbif_unmatched.py),
    in units of GBIF accepted species not already counted under the same
    'Genus epithet' by CoL in that partition.  C0 itself is unchanged."""
    if not os.path.exists(GBIF_UNMATCHED):
        die("%s is missing; run scripts/gbif_unmatched.py" % rel(GBIF_UNMATCHED))
    gb = {}
    dates = set()
    for r in read_csv(GBIF_UNMATCHED):
        gb[(r["query"], r["kingdom"])] = r
        dates.add(r["queried_on"])
    euk = [r for r in basis if is_euk(r) and as_int(r.get("c0_units"), 0) > 0]
    den = sum(as_int(r["described"], 0) for r in basis if is_euk(r))
    per = OrderedDict()
    added_all = set()
    n_q = n_res = n_missing = 0
    base_union = set()
    for r in euk:
        p = r["partition"]
        rows = shipped.col_rows.get(p, {})
        kc = Counter(x["col_kingdom"] for x in rows.values() if x["col_kingdom"])
        king = kc.most_common(1)[0][0] if kc else ""
        have = {col_canonical(x["col_name"]) for x in rows.values()
                if x["match_type"] in COL_RESOLUTION_RESOLVED and x["col_id"]}
        base_union |= {(x["col_id"]) for x in rows.values()
                       if x["match_type"] in COL_RESOLUTION_RESOLVED
                       and x["col_id"]}
        add = set()
        n_unm = 0
        for x in rows.values():
            if x["match_type"] != "unmatched" or not x["query"]:
                continue
            n_unm += 1
            g = gb.get((x["query"], king))
            n_q += 1
            if g is None:
                n_missing += 1
                continue
            if (g["matchType"] in GBIF_RESOLVED_MATCH and g["rank"] == "SPECIES"
                    and g["status"] in GBIF_RESOLVED_STATUS and g["species"]):
                n_res += 1
                if g["species"] not in have:
                    add.add(g["species"])
        if n_unm:
            per[p] = OrderedDict([
                ("unmatched_tips", n_unm),
                ("gbif_species_added", len(add)),
                ("c0_with_gbif", as_int(r["c0_units"], 0) + len(add)),
                ("percent_with_gbif",
                 pct(as_int(r["c0_units"], 0) + len(add),
                     as_int(r["described"], 0))),
            ])
        added_all |= {(p, a) for a in add}
    if n_missing:
        die("%d unmatched queries have no row in %s; run "
            "scripts/gbif_unmatched.py" % (n_missing, rel(GBIF_UNMATCHED)))
    c0 = sum(as_int(r["c0_units"], 0) for r in euk)
    num = c0 + len(added_all)
    return OrderedDict([
        ("rule", "C0 plus, per partition, the GBIF backbone accepted species "
                 "(strict match, EXACT or FUZZY, rank SPECIES) of the tips "
                 "COL26.9 leaves unmatched, not already counted under the "
                 "same 'Genus epithet'"),
        ("source", rel(GBIF_UNMATCHED)),
        ("queried_on", sorted(dates)),
        ("unmatched_tips", n_q),
        ("unmatched_tips_resolved_by_gbif", n_res),
        ("species_added", len(added_all)),
        ("numerator", num),
        ("denominator", den),
        ("percent", pct(num, den)),
        ("percent_C0", pct(c0, den)),
        ("per_partition", per),
    ])


C1B_SENSITIVITY_REASONS = ("degenerate", "at_bound", "prior_reproduced",
                           "prior_sampled")
C1B_SENSITIVITY_KINDS = ("calibration_prior_sample",)


def c1b_calibration_sensitivity(basis, shipped, truth):
    """D22, D25.  C1b excludes intervals that a calibration fixes: degenerate
    ones, ones whose centre sits on a bound, ones that reproduce the focal
    node's calibration prior (D25), and a tree set dated under root ages
    drawn from the calibration prior.  This prints C1b as it would be if they
    counted."""
    reg = crown_intervals()
    added = sorted(p for p, iv in reg.items()
                   if iv.get("c1b_exclusion") in C1B_SENSITIVITY_REASONS
                   or iv.get("kind") in C1B_SENSITIVITY_KINDS)
    alt = [dict(r) for r in basis]
    for r in alt:
        if r["partition"] in added:
            r["interval_recoverable"] = "yes"
    return OrderedDict([
        ("rule", "C1b if the intervals excluded by D22 and D25 counted"),
        ("partitions_added", added),
        ("rungs", _rung_percents(compute_rungs(alt, shipped, truth))),
    ])


# Round F3 (F3-07, F3-19).  The stricter classification that can move
# interval coverage: every counted interval whose crown node carries a
# calibration of its own (one-sided bounds included), or whose crown lies
# directly below a node with a two-sided calibration, is not counted.  A
# crown without a calibration of its own can still be held between a
# two-sided prior one node up and a minimum one node down (Amphibians:
# Data/calibration_register.csv, parent_calibration).
C1B_CROWN_CALIBRATED = ("one_sided", "two_sided", "not_checked")
C1B_PARENT_CALIBRATED = ("two_sided",)
CALIBRATION_REGISTER = os.path.join(ROOT, "Data", "calibration_register.csv")


def calibration_register():
    return {r["partition"]: r for r in read_csv(CALIBRATION_REGISTER)}


def c1b_crown_calibrated_sensitivity(basis, shipped, truth):
    reg = crown_intervals()
    cal = calibration_register()
    counted = {p: iv for p, iv in reg.items()
               if iv.get("counts_for_c1b") == "yes"}
    for p in counted:
        if (cal.get(p, {}).get("parent_calibration") or "") not in (
                "two_sided", "one_sided", "none", "not_checked"):
            die("%s: counted interval of %s has no parent_calibration"
                % (rel(CALIBRATION_REGISTER), p))
    own = sorted(p for p, iv in counted.items()
                 if iv.get("calibration") in C1B_CROWN_CALIBRATED)
    parent = sorted(p for p in counted if p not in own and cal[p][
        "parent_calibration"] in C1B_PARENT_CALIBRATED)
    drop = sorted(own + parent)
    alt = [dict(r) for r in basis]
    for r in alt:
        if r["partition"] in drop:
            r["interval_recoverable"] = "no"
    return OrderedDict([
        ("rule", "C1b with every counted interval whose crown node carries a "
                 "calibration (%s), or whose parent node carries a %s "
                 "calibration, not counted"
                 % (", ".join(C1B_CROWN_CALIBRATED),
                    " or ".join(C1B_PARENT_CALIBRATED))),
        ("partitions_removed", drop),
        ("partitions_removed_own_calibration", own),
        ("partitions_removed_parent_calibration", parent),
        ("partitions_kept", sorted(p for p in counted if p not in drop)),
        ("counted_by_crown_calibration", OrderedDict(
            (k, sorted(p for p, iv in counted.items()
                       if iv.get("calibration") == k))
            for k in ("none", "one_sided", "two_sided", "not_checked"))),
        ("rungs", _rung_percents(compute_rungs(alt, shipped, truth))),
    ])


def c1b_rule_p_sensitivity(basis, shipped, truth):
    """Round F3 (F3-14).  Rule (p) tests containment in the prior's range and
    the width ratio; it does not test where in the prior the interval lies.
    This prints C1b with the intervals that rule (p)'s test excludes counted
    (a set dated under root ages drawn from the prior, D29, stays excluded:
    its interval is the prior by construction)."""
    reg = crown_intervals()
    added = sorted(p for p, iv in reg.items()
                   if iv.get("c1b_exclusion") == "prior_reproduced"
                   and iv.get("kind") not in C1B_SENSITIVITY_KINDS)
    alt = [dict(r) for r in basis]
    for r in alt:
        if r["partition"] in added:
            r["interval_recoverable"] = "yes"
    return OrderedDict([
        ("rule", "C1b if the intervals excluded by the D25 test (p) counted"),
        ("partitions_added", added),
        ("rungs", _rung_percents(compute_rungs(alt, shipped, truth))),
    ])


def c1b_all_intervals_sensitivity(basis, shipped, truth):
    """Round-2 reader R2-06.  Interval coverage if every archived interval
    counted, whatever its kind (date assignment, taxon addition, calibration
    prior) and whether or not a calibration sets it: the upper end of the
    range that any classification of the archived intervals can give."""
    reg = crown_intervals()
    added = sorted(p for p, iv in reg.items() if iv.get("counts_for_c1b") != "yes")
    alt = [dict(r) for r in basis]
    for r in alt:
        if r["partition"] in added:
            r["interval_recoverable"] = "yes"
    return OrderedDict([
        ("rule", "C1b if every archived crown interval counted"),
        ("partitions_added", added),
        ("rungs", _rung_percents(compute_rungs(alt, shipped, truth))),
    ])


def c1_definition_sensitivity(basis, shipped, truth):
    """Audit A01.  C1 is counted as distinct accepted names among the molecular
    tips.  The superseded definition multiplied c0_units (names) by f_mol (a
    fraction of LABELS); that product assumes molecular and non-molecular
    labels resolve to, and collapse onto, accepted names at the same rate.
    Both are printed so the gap between them is always visible."""
    alt = [dict(r) for r in basis]
    for r in alt:
        r["c1_units"] = ""
    prod = _rung_percents(compute_rungs(alt, shipped, truth))
    count = _rung_percents(compute_rungs(basis, shipped, truth))
    out = OrderedDict([
        ("definition_primary",
         "C1 = distinct CoL accepted names (GTDB clusters for prokaryotes) "
         "among the tips in the molecular analysis (basis column c1_units)"),
        ("definition_product",
         "C1 = c0_units * f_mol, a label fraction applied to name units"),
        ("product", prod),
        ("difference_pp", OrderedDict(
            (k, round(count[k]["percent"] - prod[k]["percent"], PCT_DP))
            for k in ("C1", "C1d", "C1b"))),
        ("per_partition", OrderedDict(
            (r["partition"], OrderedDict([
                ("c1_count", as_float(r["c1_units"])),
                ("c1_product", round(as_int(r["c0_units"], 0)
                                     * (as_float(r["f_mol"]) or 0.0),
                                     NUMERATOR_DP)),
                ("c1_basis", r["c1_basis"]),
            ]))
            for r in basis if is_euk(r) and as_int(r.get("c0_units"), 0) > 0)),
    ])
    return out


def assembly_sensitivity(basis, shipped, truth):
    """§8.2 -- leave-assembly-class-out.

    Retaining only primary_analysis is the strict view a reader may want; it is
    printed as ONE LINE rather than being smuggled into the criterion, because
    topology_assembly is nominal and derived from unreconciled ledger fields.
    """
    classes = sorted({r["topology_assembly"] for r in basis
                      if is_euk(r) and r.get("topology_assembly")})
    out = OrderedDict()
    out["classes_present"] = OrderedDict(
        (c, sum(1 for r in basis if is_euk(r)
                and r.get("topology_assembly") == c)) for c in classes)
    for c in classes:
        keep = (lambda cc: (lambda r: r.get("topology_assembly") == cc))(c)
        out["retain_only_" + c] = OrderedDict([
            ("retained_partitions",
             sorted(r["partition"] for r in basis
                    if is_euk(r) and as_int(r.get("c0_units"), 0) > 0
                    and r.get("topology_assembly") == c)),
            ("rungs", _rung_percents(compute_rungs(basis, shipped, truth,
                                                   keep=keep))),
        ])
    for c in classes:
        keep = (lambda cc: (lambda r: r.get("topology_assembly") != cc))(c)
        out["drop_" + c] = OrderedDict([
            ("dropped_partitions",
             sorted(r["partition"] for r in basis
                    if is_euk(r) and as_int(r.get("c0_units"), 0) > 0
                    and r.get("topology_assembly") == c)),
            ("rungs", _rung_percents(compute_rungs(basis, shipped, truth,
                                                   keep=keep))),
        ])
    return out


COL_DENOMINATORS = os.path.join(ROOT, "Data", "col_denominators.csv")
PCT_REPORT_DP = 1       # D9: sensitivity percentages are reported to one decimal


def col_single_authority(basis, shipped, truth):
    """S08 / D9.  Every eukaryotic denominator replaced by the Catalogue of
    Life count of the same partition (Data/col_denominators.csv, written by
    scripts/col_denominators.py): the base release with its coverage holes
    filled from XR by the class rule stated there.  A partition with no CoL
    count at all keeps its primary denominator and is listed.  Numerators are
    unchanged; only the denominator authority changes."""
    if not os.path.exists(COL_DENOMINATORS):
        die("%s is missing; run scripts/col_denominators.py"
            % rel(COL_DENOMINATORS))
    with open(COL_DENOMINATORS, newline="", encoding="utf-8") as fh:
        col = {r["partition"]: r for r in csv.DictReader(
            ln for ln in fh if not ln.startswith("#"))}
    euk = [r for r in basis if is_euk(r)]
    missing = sorted(r["partition"] for r in euk if r["partition"] not in col)
    if missing:
        die("%s has no row for %s" % (rel(COL_DENOMINATORS), missing))
    out = OrderedDict()
    for variant, field in (("filled", "col_filled"), ("base_only", "col_base")):
        den, kept, per = 0, [], OrderedDict()
        for r in euk:
            part = r["partition"]
            c = as_int(col[part][field], 0)
            d = c if c > 0 else as_int(r["described"], 0)
            if c <= 0:
                kept.append(part)
            den += d
            n0 = as_int(r["c0_units"], 0)
            if n0:
                per[part] = OrderedDict([
                    ("described_primary", as_int(r["described"], 0)),
                    ("described_col", d),
                    ("C0_pct_primary", round(pct(n0, as_int(r["described"], 0)),
                                             PCT_REPORT_DP)),
                    ("C0_pct_col", round(pct(n0, d), PCT_REPORT_DP)),
                ])
        rungs = _rung_percents(compute_rungs(basis, shipped, truth,
                                             denominator=den))
        out[variant] = OrderedDict([
            ("denominator", den),
            ("partitions_keeping_primary", kept),
            ("rungs", rungs),
            ("rungs_pct_1dp", OrderedDict((k, round(v["percent"], PCT_REPORT_DP))
                                          for k, v in rungs.items())),
            ("per_partition", per),
        ])
    prim = _rung_percents(compute_rungs(basis, shipped, truth))
    out["primary_rungs_pct_1dp"] = OrderedDict(
        (k, round(v["percent"], PCT_REPORT_DP)) for k, v in prim.items())
    out["source"] = "%s; %s" % (rel(COL_DENOMINATORS),
                                next(iter(col.values()))["release_base"])
    return out


def denominator_sensitivity(basis, shipped, truth):
    """§8.3 -- the sensitivity the paper has never had.

    Two things are shown:
      * the represented-partitions-only denominator, which is what a reader who
        objects to counting partitions with no tree at all would use;
      * an ELASTICITY DEMONSTRATION on the single largest described-species
        count.  It perturbs that one number and shows how much every printed
        decimal moves.  IT IS NOT A CLAIM ABOUT THAT GROUP'S RICHNESS.
    """
    euk = [r for r in basis if is_euk(r)]
    den_full = sum(as_int(r["described"], 0) for r in euk)
    represented = [r for r in euk if as_int(r.get("c0_units"), 0) > 0]
    den_repr = sum(as_int(r["described"], 0) for r in represented)

    biggest = max(euk, key=lambda r: (as_int(r["described"], 0), r["partition"]))
    d_big = as_int(biggest["described"], 0)

    out = OrderedDict()
    out["full_denominator"] = OrderedDict([
        ("described", den_full),
        ("n_partitions", len(euk)),
        ("rungs", _rung_percents(compute_rungs(basis, shipped, truth))),
    ])
    out["represented_partitions_only"] = OrderedDict([
        ("described", den_repr),
        ("n_partitions", len(represented)),
        ("note", "the denominator restricted to partitions that have a "
                 "canonical tree at all; the structural zeros are removed from "
                 "BOTH sides"),
        ("rungs", _rung_percents(compute_rungs(basis, shipped, truth,
                                               denominator=den_repr))),
    ])
    out["col_single_authority"] = col_single_authority(basis, shipped, truth)
    out["largest_single_count"] = OrderedDict([
        ("partition", biggest["partition"]),
        ("described", d_big),
        ("described_source", biggest.get("described_source")),
        ("share_of_denominator_pct", pct(d_big, den_full)),
        ("note", "one round number from a single compilation carries this "
                 "share of the whole eukaryotic denominator; every printed "
                 "decimal in the paper is elastic to it"),
    ])
    elasticity = OrderedDict()
    for sign, tag in ((-1.0, "minus"), (1.0, "plus")):
        delta = sign * ELASTICITY_FRACTION * d_big
        den_alt = den_full + delta
        elasticity["%s_%d_pct" % (tag, int(round(100 * ELASTICITY_FRACTION)))] = \
            OrderedDict([
                ("denominator", int(round(den_alt))),
                ("rungs", _rung_percents(compute_rungs(
                    basis, shipped, truth, denominator=den_alt))),
            ])
    out["elasticity_demonstration"] = OrderedDict([
        ("perturbed_partition", biggest["partition"]),
        ("fraction", ELASTICITY_FRACTION),
        ("warning", "AN ELASTICITY DEMONSTRATION, NOT A RICHNESS CLAIM. The "
                    "perturbation is a device for showing how much of every "
                    "printed decimal rests on one round number."),
        ("results", elasticity),
    ])
    return out


# ==========================================================================
# DISCLOSURES the Results MUST print
# ==========================================================================

def disclosures(basis, shipped, truth, mol, rungs, warnings):
    euk = [r for r in basis if is_euk(r)]
    den = sum(as_int(r["described"], 0) for r in euk)
    represented = [r for r in euk if as_int(r.get("c0_units"), 0) > 0]

    # -- dominance of the single largest contributor -----------------------
    def contrib(row, rung):
        return as_int(row["c0_units"], 0) * rung_weight(row, rung)

    top = max(represented, key=lambda r: (contrib(r, "C0"), r["partition"]))
    dominance = OrderedDict([("partition", top["partition"])])
    for rung in RUNGS:
        tot = rungs[rung]["numerator_sum"]
        # Each share ships the EXACT pair that produces it: this partition's
        # (fractional) contribution over the rung total.  Without the pair the
        # share is an unreachable percentage.
        dominance["contribution_to_%s" % rung] = round(contrib(top, rung),
                                                       NUMERATOR_DP)
        dominance["%s_total" % rung] = tot
        dominance["share_of_%s_pct" % rung] = pct(contrib(top, rung), tot)
    dominance["note"] = ("state this and frame the asymmetry as the finding, "
                         "not as a caveat")

    # -- f_basis accounting of C1 -------------------------------------------
    by_basis = defaultdict(float)
    parts_by_basis = defaultdict(list)
    for r in represented:
        by_basis[r["f_basis"]] += contrib(r, "C1")
        parts_by_basis[r["f_basis"]].append(r["partition"])
    c1 = rungs["C1"]["numerator_sum"]
    accounting = OrderedDict()
    for b in sorted(by_basis, key=lambda k: (-by_basis[k], k)):
        accounting[b] = OrderedDict([
            ("tip_equivalents", int(round(by_basis[b]))),
            ("tip_equivalents_exact", round(by_basis[b], NUMERATOR_DP)),
            ("C1_total", c1),
            ("share_of_C1_pct", pct(by_basis[b], c1)),
            ("partitions", sorted(parts_by_basis[b])),
        ])

    # -- overlapping labels --------------------------------------------------
    dupes, _ = shared_label_detail(basis, shipped)
    overlap = OrderedDict([
        ("n_labels_in_two_or_more_partitions", len(dupes)),
        ("share_of_C0_union_pct", pct(len(dupes), rungs["C0"]["numerator_union"])),
        ("examples", OrderedDict(
            (lab, dupes[lab]) for lab in sorted(dupes))),
        ("note", "outgroups and homonyms. C0 is quoted on the UNION "
                 "convention; the per-partition SUM is given alongside so the "
                 "two can never drift. EVERY derived quantity -- dark matter "
                 "above all -- must use the union figure consistently."),
    ])

    # -- structural zeros -----------------------------------------------------
    zeros = [r for r in euk if as_int(r.get("c0_units"), 0) == 0]
    structural = OrderedDict([
        ("n_partitions", len(zeros)),
        ("of_eukaryotic_partitions", len(euk)),
        ("described", sum(as_int(r["described"], 0) for r in zeros)),
        ("share_of_denominator_pct",
         pct(sum(as_int(r["described"], 0) for r in zeros), den)),
        ("partitions", sorted(r["partition"] for r in zeros)),
        ("note", "a STRUCTURAL zero, not a criterion-dependent one -- and now "
                 "computable from the deposit file itself, because these "
                 "partitions are rows of coverage_basis.csv."),
    ])

    # -- dark matter, on the union convention throughout ---------------------
    c0u = rungs["C0"]["numerator_union"]
    dark = OrderedDict([
        ("no_tip_anywhere", den - c0u),
        ("no_tip_anywhere_pct", pct(den - c0u, den)),
        ("no_molecular_placement", den - int(round(c1))),
        ("no_molecular_placement_exact", round(den - c1, NUMERATOR_DP)),
        ("no_molecular_placement_pct", pct(den - c1, den)),
        ("convention", "formed on the C0 UNION numerator (%d), never on the "
                       "per-partition sum" % c0u),
        ("replaces", "the retired 23.7 % figure, which was a hard-coded "
                     "literal reachable from no numerator/denominator pair in "
                     "the deposit. WITHDRAWN, not adjusted."),
    ])

    # -- the ordering fix: what the old strict criterion could not see -------
    ordering = []
    for r in sorted(represented,
                    key=lambda x: (-(as_float(x.get("f_mol")) or 0.0),
                                   x["partition"])):
        d = as_int(r["described"], 0)
        # Each row carries its own numerator and denominator, so every
        # percentage in the table that justifies the redesign is reachable
        # from the row that prints it -- no reader has to reconstruct a
        # denominator, which is how 23.7 % survived three submissions.
        ordering.append(OrderedDict([
            ("partition", r["partition"]),
            ("f_mol", as_float(r.get("f_mol"))),
            ("topology_assembly", r.get("topology_assembly")),
            ("described", d),
            ("denominator", d),
            ("c0_units", as_int(r["c0_units"], 0)),
            ("C0_pct", pct(as_int(r["c0_units"], 0), d)),
            ("C1_tip_equivalents", round(contrib(r, "C1"), NUMERATOR_DP)),
            ("C1_pct", pct(contrib(r, "C1"), d)),
            ("old_strict_pct",
             0.0 if r.get("topology_assembly") != ASSEMBLY_PRIMARY
             else pct(contrib(r, "C1"), d)),
        ]))

    return OrderedDict([
        ("dominance", dominance),
        ("f_basis_accounting_of_C1", accounting),
        ("overlapping_labels", overlap),
        ("structural_zeros", structural),
        ("dark_matter", dark),
        ("ordering_fix", ordering),
    ])


# ==========================================================================
# BUILDERS
# ==========================================================================

def build_coverage_basis(basis, warnings):
    return render_csv(BASIS_HEADER, basis)


def build_molecular_tips(truth, mol, shipped, warnings):
    """(partition, standardized_id) for every partition whose molecular source
    is a deposited tree.  For those partitions C1 is a genuine enumerable,
    dedupable, prunable TIP SET rather than an expectation.

    Until EVERY partition is enumerated here, the paper writes
    "tip-equivalents", never "species", at C1 and below."""
    rows = []
    for part in sorted(mol):
        ev = mol[part]
        if ev.tip_names is None or not ev.source_file:
            continue
        idmap = shipped.ids.get(part, {})
        for nm in sorted(ev.tip_names):
            sid = idmap.get(nm)
            if sid is None:
                continue
            rows.append(OrderedDict([("partition", part),
                                     ("standardized_id", sid)]))
    if not rows:
        warnings.append(
            "molecular_tips.csv is EMPTY: no partition resolved to a deposited "
            "molecular source tree. C1 is then an expectation everywhere and "
            "the paper must say 'tip-equivalents' throughout.")
    return render_csv(["partition", "standardized_id"], rows)


def build_partition_summary(truth, pt, warnings):
    """Data/site/partition_summary.csv -- the web app's per-partition ledger.
    Column set preserved exactly as deposited so the site build keeps working.
    `dated` follows the SHIPPED FILE's geometry, never a ledger flag."""
    header = ["partition", "group", "category", "phylo_number",
              "described_species", "estimated_total", "estimated_low",
              "estimated_high", "is_unrepresented", "dated", "study",
              "filename"]
    rows = []
    for part, sp in truth.species.items():
        unrep = truth.is_unrepresented(part)
        rows.append(OrderedDict([
            ("partition", part),
            ("group", (sp.get("group") or "").strip()),
            ("category", (sp.get("category") or "").strip()),
            ("phylo_number", "" if unrep else truth.labels_unique(part)),
            ("described_species", truth.described(part)),
            ("estimated_total", pt.estimated_total(part)),
            ("estimated_low", as_int(sp.get("estimated_low"), "")),
            ("estimated_high", as_int(sp.get("estimated_high"), "")),
            ("is_unrepresented", "True" if unrep else "False"),
            ("dated", "False" if unrep
             else ("True" if truth.is_chronogram(part) else "False")),
            ("study", "" if unrep
             else (truth.canon_prov[part].get("citation_key") or "").strip()),
            ("filename", "" if unrep else truth.canon_file[part]),
        ]))
    return render_csv(header, rows)


def build_coverage_per_partition(basis, truth, warnings):
    """Paper/inputs/coverage_per_partition.csv -- THE per-partition ladder table.

    THE DATING FIX.  Dating is a FILTER, not a coverage column here.  The
    columns are C0 %, C1 %, is_chronogram, roottotip_cv, interval_kind, so a
    completely sampled but uncalibrated partition reads
    "91.39 / 91.39 / not ultrametric (CV 0.296) / --" -- the true statement --
    instead of "91.39 / 91.39 / 0 / 0", which is not.  The nested chain
    survives in the abstract sentence, where it is a legitimate refinement of
    one proposition.
    """
    header = ["partition", "category", "described", "c0_units",
              "f_mol", "f_low", "f_high", "f_basis",
              "coverage_C0_pct", "coverage_C1_pct",
              "is_chronogram", "roottotip_cv", "interval_kind", "posterior_n",
              "topology_assembly"]
    rows = []
    for r in basis:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        d = as_int(r["described"], 0)
        tips = as_int(r["c0_units"], 0)
        if not d:
            warnings.append("partition %r has described == 0; its coverage "
                            "percentages are emitted blank" % r["partition"])
        rows.append(OrderedDict([
            ("partition", r["partition"]),
            ("category", r["category"]),
            ("described", d),
            ("c0_units", tips),
            ("f_mol", r["f_mol"]),
            ("f_low", r["f_low"]),
            ("f_high", r["f_high"]),
            ("f_basis", r["f_basis"]),
            ("coverage_C0_pct", fnum(pct(tips, d), PCT_DP)),
            ("coverage_C1_pct",
             fnum(pct(tips * rung_weight(r, "C1"), d), PCT_DP)),
            ("is_chronogram", r["is_chronogram"]),
            ("roottotip_cv", r["roottotip_cv"]),
            ("interval_kind", r["interval_kind"] or ""),
            ("posterior_n", r["posterior_n"]),
            ("topology_assembly", r["topology_assembly"]),
        ]))
    rows.sort(key=lambda x: (-(as_float(x["coverage_C0_pct"]) or 0.0),
                             x["partition"]))
    return render_csv(header, rows)


def build_fig1_coverage_bars(truth, pt, warnings):
    """Figure 1 -- the ALL-LIFE inventory view (all partitions, prokaryotes
    included).  This is the one place the all-life described total is
    legitimately used, and the caption must say so.  in_dated_tree /
    in_undated_tree route the adopted tip count by the SHIPPED FILE's
    ultrametricity, never by a ledger flag."""
    header = ["group", "category", "estimated_total", "estimated_low",
              "estimated_high", "described", "in_dated_tree", "in_undated_tree",
              "estimate_source", "confidence", "is_unrepresented", "study",
              "tips", "dated"]
    order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    ordered = sorted(truth.species.items(),
                     key=lambda kv: (order.get(kv[1].get("category"), len(order)),
                                     -(truth.described(kv[0]) or 0), kv[0]))
    rows = []
    for part, sp in ordered:
        unrep = truth.is_unrepresented(part)
        tips = "" if unrep else truth.labels_unique(part)
        chrono = (not unrep) and truth.is_chronogram(part)
        rows.append(OrderedDict([
            ("group", part),
            ("category", (sp.get("category") or "").strip()),
            ("estimated_total", pt.estimated_total(part)),
            ("estimated_low", as_int(sp.get("estimated_low"), "")),
            ("estimated_high", as_int(sp.get("estimated_high"), "")),
            ("described", truth.described(part)),
            ("in_dated_tree", 0 if unrep else (tips if chrono else 0)),
            ("in_undated_tree", 0 if unrep else (0 if chrono else tips)),
            ("estimate_source", (sp.get("estimate_source") or "").strip()),
            ("confidence", pt.confidence(part)),
            ("is_unrepresented", "True" if unrep else "False"),
            ("study", "" if unrep
             else (truth.canon_prov[part].get("citation_key") or "").strip()),
            ("tips", tips),
            ("dated", "False" if unrep else ("True" if chrono else "False")),
        ]))
    return render_csv(header, rows)


def build_fig2_scatter(basis, warnings):
    """Figure 2 -- REPLACES fig2_sensitivity_bounds.py entirely.

    Per-partition scatter: x = f_mol, y = C0 coverage %, marker shape =
    topology_assembly, marker area proportional to described species.  This is
    the only artefact that makes the assembly fact VISIBLE rather than
    asserted: a large, highly covered, automated-synthesis partition plots at
    the top right in its own shape, and a huge, barely covered one plots at the
    bottom left.
    """
    header = ["partition", "category", "described", "c0_units",
              "f_mol", "f_low", "f_high", "f_basis",
              "coverage_C0_pct", "coverage_C1_pct", "topology_assembly",
              "is_chronogram", "label"]
    rows = []
    for r in basis:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        d = as_int(r["described"], 0)
        tips = as_int(r["c0_units"], 0)
        rows.append(OrderedDict([
            ("partition", r["partition"]),
            ("category", r["category"]),
            ("described", d),
            ("c0_units", tips),
            ("f_mol", r["f_mol"]),
            ("f_low", r["f_low"]),
            ("f_high", r["f_high"]),
            ("f_basis", r["f_basis"]),
            ("coverage_C0_pct", fnum(pct(tips, d), PCT_DP)),
            ("coverage_C1_pct",
             fnum(pct(tips * rung_weight(r, "C1"), d), PCT_DP)),
            ("topology_assembly", r["topology_assembly"]),
            ("is_chronogram", r["is_chronogram"]),
            ("label", r["partition"]),
        ]))
    rows.sort(key=lambda x: (-x["described"], x["partition"]))
    return render_csv(header, rows)


def build_fig3_c0c1_bars(basis, warnings):
    """Figure 3 -- two panels, C0 and C1 per partition, sorted, WITH the
    treeless partitions shown as STRUCTURAL ZEROS in a distinct hatch.  They
    are rows here (is_structural_zero = yes) precisely so the figure cannot
    quietly omit them."""
    header = ["partition", "category", "described", "c0_units",
              "coverage_C0_pct", "coverage_C1_pct", "is_structural_zero",
              "is_chronogram", "topology_assembly"]
    rows = []
    for r in basis:
        if not is_euk(r):
            continue
        d = as_int(r["described"], 0)
        tips = as_int(r["c0_units"], 0)
        zero = tips == 0
        rows.append(OrderedDict([
            ("partition", r["partition"]),
            ("category", r["category"]),
            ("described", d),
            ("c0_units", tips),
            ("coverage_C0_pct", fnum(pct(tips, d), PCT_DP)),
            ("coverage_C1_pct",
             fnum(pct(tips * rung_weight(r, "C1"), d), PCT_DP)),
            ("is_structural_zero", "yes" if zero else "no"),
            ("is_chronogram", r["is_chronogram"] or ""),
            ("topology_assembly", r["topology_assembly"] or ""),
        ]))
    rows.sort(key=lambda x: (-(as_float(x["coverage_C0_pct"]) or 0.0),
                             x["partition"]))
    return render_csv(header, rows)


def build_fig3_counts(truth, warnings):
    """Counts of the partition-canonical trees by category and by the SHIPPED
    FILE's ultrametricity verdict.  Under the settled rule this moves several
    partitions out of "Dated": the flag now follows the file."""
    header = ["category", "type", "count"]
    tally = defaultdict(lambda: {"Dated": 0, "Undated": 0})
    for part in truth.canon_prov:
        cat = (truth.species[part].get("category") or "").strip()
        tally[cat]["Dated" if truth.is_chronogram(part) else "Undated"] += 1
    cats = ([c for c in CATEGORY_ORDER if c in tally]
            + sorted(c for c in tally if c not in CATEGORY_ORDER))
    rows = []
    for cat in cats:
        for kind in ("Dated", "Undated"):
            rows.append(OrderedDict([("category", cat), ("type", kind),
                                     ("count", tally[cat][kind])]))
    return render_csv(header, rows)


# -------------------------------------------------------------- Table S5 ----

S5_HEADER = ["partition_name", "study", "doi", "tips_unique", "tips_raw",
             "is_chronogram", "roottotip_cv", "tree_class", "f_mol", "f_basis",
             "molecular_source_file", "topology_assembly", "placement_evidence",
             "imputation_status", "dating_code", "archival_class",
             "interval_kind", "posterior_n", "source_of_truth", "notes",
             "conflict_note"]


def _s5_rows(truth, basis, mol, ledger):
    """Supplementary Table S5 -- canonical tree provenance, rebuilt on the
    measured columns.  `dated` is gone; is_chronogram + roottotip_cv replace it,
    so a reader can apply their own ultrametricity tolerance."""
    by_part = {r["partition"]: r for r in basis}
    rows = []
    for part in sorted(truth.canon_prov):
        pv = truth.canon_prov[part]
        b = by_part[part]
        rows.append(OrderedDict([
            ("partition_name", part),
            ("study", (pv.get("citation_key") or "").strip()),
            ("doi", (pv.get("doi") or "").strip()),
            ("tips_unique", truth.labels_unique(part)),
            ("tips_raw", truth.labels_raw(part)),
            ("is_chronogram", b["is_chronogram"]),
            ("roottotip_cv", b["roottotip_cv"]),
            ("tree_class", ledger.tree_class(part)),
            ("f_mol", b["f_mol"]),
            ("f_basis", b["f_basis"]),
            ("molecular_source_file", b["molecular_source_file"]),
            ("topology_assembly", b["topology_assembly"]),
            ("placement_evidence", b["placement_evidence"]),
            ("imputation_status", (pv.get("imputation_status") or "").strip()),
            ("dating_code", (pv.get("dating_method_code") or "").strip()),
            ("archival_class", (pv.get("archival_class") or "").strip()),
            ("interval_kind", b["interval_kind"]),
            ("posterior_n", b["posterior_n"]),
            ("source_of_truth", (pv.get("source_of_truth") or "").strip()),
            ("notes", ledger.notes(part)),
            ("conflict_note", (pv.get("conflict_note") or "").strip()),
        ]))

    # The aggregate family-level block that truth_provenance carries as a
    # single non-canonical row, with its dated count RECOMPUTED over the actual
    # shipped files, minus any file that is itself a partition canonical.
    canon_files = set(truth.canon_file.values())
    agg = next((r for r in truth.alt_prov
                if "*" in (r.get("filename") or "")), None)
    if agg is not None:
        stem = (agg.get("filename") or "").split("*")[0]
        block = sorted(f for f in truth.dated
                       if stem and f.startswith(stem) and f not in canon_files)
        if block:
            n_dated = sum(1 for f in block
                          if truthy(truth.dated[f]["is_chronogram"]))
            named_total = len([f for f in truth.dated
                               if stem and f.startswith(stem)])
            rows.append(OrderedDict([
                # The label is the truth file's own; only the COUNT of trees
                # actually aggregated here is recomputed, because any file in
                # the block that is also a partition canonical is reported in
                # its own row and must not be counted twice.
                ("partition_name", (agg.get("partition_name") or "").strip()
                 or "aggregate family-level block"),
                ("study", (agg.get("citation_key") or "").strip()),
                ("doi", (agg.get("doi") or "").strip()),
                ("tips_unique", sum(as_int(truth.tips[f]["tips_adopted"], 0)
                                    for f in block)),
                ("tips_raw", sum(as_int(truth.tips[f]["raw_leaves"], 0)
                                 for f in block)),
                ("is_chronogram", "%d of %d" % (n_dated, len(block))),
                ("roottotip_cv", ""),
                ("tree_class", "family level compilation"),
                ("f_mol", fnum(as_float(agg.get("molecular_tip_fraction"), 1.0),
                               FRAC_DP)),
                ("f_basis", "aggregate_row"),
                ("molecular_source_file", ""),
                ("topology_assembly", ""),
                ("placement_evidence", ""),
                ("imputation_status", (agg.get("imputation_status") or "").strip()),
                ("dating_code", (agg.get("dating_method_code") or "").strip()),
                ("archival_class", (agg.get("archival_class") or "").strip()),
                ("interval_kind", ""),
                ("posterior_n", (agg.get("posterior_n") or "").strip()),
                ("source_of_truth", (agg.get("source_of_truth") or "").strip()),
                ("notes", "ONE row aggregates %d of the %d family-level trees "
                          "in this block. The %d file(s) that are also "
                          "partition canonicals are reported in their own rows "
                          "above and are excluded here to avoid "
                          "double-counting; the label's own count is the "
                          "block total, not the aggregated count."
                          % (len(block), named_total, named_total - len(block))),
                ("conflict_note", (agg.get("conflict_note") or "").strip()),
            ]))
    return rows


def build_table_s5_csv(truth, basis, mol, ledger, warnings):
    return render_csv(S5_HEADER, _s5_rows(truth, basis, mol, ledger))


def build_table_s5_tex(truth, basis, mol, ledger, warnings):
    cols = (r"@{}p{1.5cm}p{1.8cm}r c c p{0.9cm} p{1.6cm} p{1.5cm} "
            r"p{0.9cm} p{3.0cm}@{}")
    head = (r"\textbf{Partition} & \textbf{Study} & \textbf{Tips} & "
            r"\textbf{Chrono.} & \textbf{CV} & \textbf{$f_{\mathrm{mol}}$} & "
            r"\textbf{$f$ basis} & \textbf{Assembly} & \textbf{Dating} & "
            r"\textbf{Notes} \\")
    body = []
    for r in _s5_rows(truth, basis, mol, ledger):
        cv = r["roottotip_cv"]
        body.append(" & ".join([
            tex_escape(r["partition_name"]),
            tex_escape(r["study"]),
            thousands(r["tips_unique"]),
            tex_escape(r["is_chronogram"]),
            tex_escape(("%.4f" % as_float(cv)) if as_float(cv) is not None else "--"),
            tex_escape(("%.3f" % as_float(r["f_mol"]))
                       if as_float(r["f_mol"]) is not None else "--"),
            tex_escape(r["f_basis"].replace("_", " ")),
            tex_escape((r["topology_assembly"] or "--").replace("_", " ")),
            tex_escape(r["dating_code"]),
            tex_escape(r["notes"]),
        ]) + r" \\")
    return _longtable(cols, head, body, small=True)


def _longtable(colspec, header_row, body_lines, small=False, caption=None):
    """A self-contained longtable fragment for \\input."""
    parts = []
    if small:
        parts += [r"\begingroup", r"\scriptsize",
                  r"\setlength{\tabcolsep}{3pt}",
                  r"\setlength{\LTleft}{0pt}", r"\setlength{\LTright}{0pt}"]
    parts.append(r"\begin{longtable}{%s}" % colspec)
    if caption:
        parts.append(r"\caption{%s} \\" % caption)
    parts += [r"\toprule", header_row, r"\midrule", r"\endfirsthead",
              r"\toprule", header_row, r"\midrule", r"\endhead",
              r"\midrule", r"\endfoot",
              r"\bottomrule", r"\endlastfoot"]
    parts += body_lines
    parts.append(r"\end{longtable}")
    if small:
        parts.append(r"\endgroup")
    return "\n".join(parts) + "\n"


# ------------------------------------------------- Table S6 (two readings) ---
#
# The deposit uses the label "S6" for two different tables:
#   supplementary_body.tex -> "S6: Source data for Figure 1"
#   manuscript_body.tex    -> cites an alternative-trees table
# Both are regenerated under distinct filenames and the collision is FLAGGED.

S6_FIG1_HEADER = ["group", "category", "estimated_total", "described",
                  "in_dated_tree", "in_undated_tree", "source"]


def _s6_fig1_rows(truth, pt):
    order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    ordered = sorted(
        truth.species.items(),
        key=lambda kv: (truth.is_unrepresented(kv[0]),
                        order.get(kv[1].get("category"), len(order)),
                        -(truth.described(kv[0]) or 0), kv[0]))
    rows = []
    for part, sp in ordered:
        unrep = truth.is_unrepresented(part)
        tips = None if unrep else truth.labels_unique(part)
        chrono = (not unrep) and truth.is_chronogram(part)
        rows.append(OrderedDict([
            ("group", part),
            ("category", "Not yet represented" if unrep
             else (sp.get("category") or "").strip()),
            ("estimated_total", pt.estimated_total(part)),
            ("described", truth.described(part)),
            ("in_dated_tree", tips if chrono else ""),
            ("in_undated_tree", "" if (unrep or chrono) else tips),
            ("source", (sp.get("estimate_source") or "").strip()),
        ]))
    return rows


def build_table_s6_fig1_csv(truth, pt, warnings):
    return render_csv(S6_FIG1_HEADER, _s6_fig1_rows(truth, pt))


def build_table_s6_fig1_tex(truth, pt, warnings):
    cols = r"p{2.2cm}p{2.5cm}rrrrl"
    head = (r"\textbf{Group} & \textbf{Category} & \textbf{Est.\ total} & "
            r"\textbf{Described} & \textbf{Dated} & \textbf{Undated} & "
            r"\textbf{Source} \\")
    body, prev = [], None
    for r in _s6_fig1_rows(truth, pt):
        if prev is not None and r["category"] != prev:
            body.append(r"\midrule")
        prev = r["category"]
        body.append(" & ".join([
            tex_escape(r["group"]), tex_escape(r["category"]),
            thousands(r["estimated_total"]), thousands(r["described"]),
            thousands(r["in_dated_tree"]) if r["in_dated_tree"] != "" else "--",
            thousands(r["in_undated_tree"]) if r["in_undated_tree"] != "" else "--",
            tex_escape(r["source"]),
        ]) + r" \\")
    return _longtable(cols, head, body)



# ------------------------------------------------ the dataset table (Table) ---
#
# F01-F07: the dataset table of the main text is GENERATED, so its tip counts,
# accepted-name counts, dating column and list of unrepresented partitions are
# the basis table's own values and cannot disagree with it.  The only curated
# input is Data/dataset_table_registry.csv: per shipped file, the display
# label, BibTeX key, publication year and journal (literature metadata), and,
# for NON-canonical trees only, a described-species count with its source.
# Canonical rows take `described` and `described_source` from the basis table.
DATASET_REGISTRY = os.path.join(ROOT, "Data", "dataset_table_registry.csv")
# Round F3 (F3-06): larger trees of a group or a subclade that the May 2026
# search missed and that are not counted (decision D40); Table 1 marks each
# group with one and its caption names the tree and its size as the source
# states it.
LARGER_TREES = os.path.join(ROOT, "Data", "larger_trees_not_counted.csv")


def larger_trees_not_counted():
    out = defaultdict(list)
    for r in read_csv(LARGER_TREES):
        out[r["partition"]].append(r)
    return out


def larger_tree_note():
    """The Table 1 caption note: per group, each larger tree not counted."""
    parts = []
    for p, rs in sorted(larger_trees_not_counted().items()):
        items = []
        for r in rs:
            n = ("%s " % r["n_qualifier"] if r["n_qualifier"] else "") \
                + thousands(r["n"])
            items.append(r"\citet{%s}, %s %s of %s%s"
                         % (r["bibkey"], n, r["unit"], r["clade"],
                            ", dated" if r["dated"] == "yes" else ""))
        parts.append("%s: %s" % (p, " and ".join(items)))
    return "; ".join(parts)


REFERENCES_BIB = os.path.join(ROOT, "Paper", "references.bib")
DATASET_SECTION_CROSS = "Cross-cutting"
DATASET_HEADER = ["section", "group", "filename", "bibkey", "year", "journal",
                  "is_canonical", "tips", "accepted_names", "described",
                  "coverage_pct", "dated", "roottotip_cv", "described_source",
                  "flags"]


def _bib_years():
    years = {}
    if not os.path.exists(REFERENCES_BIB):
        return years
    text = open(REFERENCES_BIB, encoding="utf-8").read()
    for m in re.finditer(r"^@\w+\{([^,]+),(.*?)^\}", text, re.M | re.S):
        y = re.search(r"\byear\s*=\s*[{\"]?(\d{4})", m.group(2))
        if y:
            years[m.group(1).strip()] = int(y.group(1))
    return years


def _noncanonical_accepted_names():
    out = defaultdict(set)
    for r in read_csv(COL_RESOLUTION_NONCANON, required=False):
        out[r["filename"]]          # a file with rows but no match counts 0
        if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
            out[r["filename"]].add(r["col_id"])
    return out


def _dataset_rows(truth, basis, warnings):
    reg = read_csv(DATASET_REGISTRY)
    years = _bib_years()
    by_tree = {r["canonical_tree"]: r for r in basis if r.get("canonical_tree")}
    by_part = {r["partition"]: r for r in basis}
    nonc = _noncanonical_accepted_names()
    scales = time_scales()
    cat_order = ["Vertebrates", "Plants", "Arthropods", "Other animals",
                 "Microbes & protists", DATASET_SECTION_CROSS]

    listed = set()
    rows = []
    for g in reg:
        fn = g["filename"].strip()
        flags = set(f for f in (g.get("flags") or "").split(";") if f)
        if g["bibkey"] not in years:
            die("dataset registry: %s cites %s, which is not in %s"
                % (fn, g["bibkey"], rel(REFERENCES_BIB)))
        if as_int(g["year"], 0) != years[g["bibkey"]]:
            die("dataset registry: %s gives year %s but %s says %d for %s"
                % (fn, g["year"], rel(REFERENCES_BIB), years[g["bibkey"]],
                   g["bibkey"]))
        row = OrderedDict((k, "") for k in DATASET_HEADER)
        row.update({"group": g["group_label"], "filename": fn,
                    "bibkey": g["bibkey"], "year": g["year"],
                    "journal": g["journal"], "flags": ";".join(sorted(flags))})
        if "aggregate" in flags:
            pat = re.compile("^" + re.escape(fn).replace(r"\*", ".*") + "$")
            files = [f for f in truth.dated if pat.match(f)]
            if not files:
                die("dataset registry: pattern %s matches no shipped file" % fn)
            listed.update(files)
            n_chrono = sum(truthy(truth.dated[f]["is_chronogram"])
                           for f in files)
            row.update({"section": DATASET_SECTION_CROSS, "is_canonical": "no",
                        "tips": "%d trees" % len(files),
                        "dated": "%d of %d" % (n_chrono, len(files))})
            rows.append(row)
            continue
        if "cited_not_distributed" in flags:
            row.update({"section": DATASET_SECTION_CROSS, "is_canonical": "no"})
            rows.append(row)
            continue
        if fn not in truth.dated:
            die("dataset registry lists %s, which is not a shipped tree" % fn)
        listed.add(fn)
        t = truth.dated[fn]
        chrono = truthy(t["is_chronogram"])
        b = by_tree.get(fn)
        if b is not None:
            if (g.get("described") or "").strip():
                die("dataset registry: %s is canonical; its described count "
                    "comes from the basis table, not the registry" % fn)
            section = b["category"]
            names = as_int(b["c0_units"], 0)
            desc = as_int(b["described"], 0)
            src = b["described_source"]
            if truthy(b.get("is_prokaryotic")):
                flags.add("gtdb_clusters")
            scale = scales.get(fn, "")
            dated = ("yes" if b["is_chronogram"] == "yes"
                     else "relative" if scale == "relative" else "no")
            canon = "yes"
        else:
            part = (t.get("partition_name") or "").strip()
            section = (by_part[part]["category"] if part in by_part
                       else DATASET_SECTION_CROSS)
            names = len(nonc.get(fn, ())) if fn in nonc else ""
            desc = as_int(g.get("described"), 0) or ""
            src = g.get("described_source") or ""
            if desc and not src:
                die("dataset registry: %s has a described count without a "
                    "source" % fn)
            dated = "ultrametric" if chrono else "no"
            canon = "no"
        cov = (fnum(100.0 * names / desc, 1)
               if desc and names != "" else "")
        row.update({"section": section, "is_canonical": canon,
                    "tips": as_int(t["distinct_labels"], 0),
                    "accepted_names": names, "described": desc,
                    "coverage_pct": cov, "dated": dated,
                    "roottotip_cv": t["roottotip_cv"],
                    "described_source": src,
                    "flags": ";".join(sorted(flags))})
        rows.append(row)

    missing = sorted(f for f in truth.dated if f not in listed)
    if missing:
        die("dataset registry does not list %d shipped tree(s): %s"
            % (len(missing), ", ".join(missing[:8])))
    for r in basis:
        if is_euk(r) and r.get("is_unrepresented") == "yes":
            row = OrderedDict((k, "") for k in DATASET_HEADER)
            row.update({"section": "Not yet represented",
                        "group": r["partition"], "is_canonical": "no",
                        "tips": 0, "accepted_names": 0,
                        "described": as_int(r["described"], 0),
                        "coverage_pct": fnum(0.0, 1),
                        "described_source": r["described_source"]})
            rows.append(row)
    cat_order.append("Not yet represented")

    def key(r):
        sec = cat_order.index(r["section"]) if r["section"] in cat_order \
            else len(cat_order)
        canon = 0 if r["is_canonical"] == "yes" else 1
        size = (-as_int(r["described"], 0) if r["section"] ==
                "Not yet represented" else
                -as_int(r["tips"], 0) if isinstance(r["tips"], int) else 0)
        return (sec, canon, size, r["group"])
    rows.sort(key=key)
    return rows


def build_dataset_table_csv(truth, basis, warnings):
    return render_csv(DATASET_HEADER, _dataset_rows(truth, basis, warnings))


DATASET_DATED_TEX = {"yes": "Yes", "no": "No", "relative": "Rel.$^f$",
                     "ultrametric": "Ultr.$^g$"}


def build_dataset_table_tex(truth, basis, warnings):
    """Body rows for the main-text dataset table (10 columns); the longtable
    environment, header and caption live in the wrappers."""
    out, prev = [], None
    for r in _dataset_rows(truth, basis, warnings):
        if r["section"] != prev:
            if prev is not None:
                out.append(r"\midrule")
            out.append(r"\textbf{\textit{%s}} & & & & & & & & & \\"
                       % tex_escape(r["section"]))
            prev = r["section"]
        flags = set(r["flags"].split(";")) if r["flags"] else set()
        group = tex_escape(r["group"]).replace(r"\textbackslash{}", "\\")
        if r["is_canonical"] == "yes":
            group = r"\textbf{%s}" % group
        if r["bibkey"]:
            study = r"\citet{%s}" % r["bibkey"]
            journal = tex_escape(r["journal"])
            if "preprint" in flags:
                journal += "$^e$"
            year = r["year"]
        else:
            study = journal = year = "---"
        tips = r["tips"]
        if isinstance(tips, int) or (isinstance(tips, str) and tips.isdigit()):
            tips = thousands(tips)
            if "genus_level" in flags:
                tips += "$^b$"
        elif tips == "":
            tips = "---"
        names = (thousands(r["accepted_names"]) if r["accepted_names"] != ""
                 else "---")
        desc = thousands(r["described"]) if r["described"] != "" else "---"
        if "gtdb_clusters" in flags:
            names += "$^c$"
            desc += "$^c$"
        cov = r["coverage_pct"] if r["coverage_pct"] != "" else "---"
        if "cited_not_distributed" in flags:
            tips = names = desc = cov = "---"
            dated = "---$^d$"
        elif "aggregate" in flags:
            dated = tex_escape(r["dated"])
        else:
            dated = DATASET_DATED_TEX.get(r["dated"], "---")
        src = tex_escape(r["described_source"]) or "---"
        out.append(" & ".join([group, study, str(year), journal, tips, names,
                               desc, cov, dated, src]) + r" \\")
    return "\n".join(out) + "\n"


def build_dataset_table_macros(truth, basis, warnings):
    rows = _dataset_rows(truth, basis, warnings)
    canon = [r for r in rows if r["is_canonical"] == "yes"]
    unrep = [r for r in rows if r["section"] == "Not yet represented"]
    shipped = [r for r in rows if r["filename"] in truth.dated]
    agg = [r for r in rows if "aggregate" in r["flags"]]
    n_agg = sum(int(r["tips"].split()[0]) for r in agg)
    euk_d = sum(as_int(r["described"], 0) for r in basis if is_euk(r))
    lines = [
        "% generated by scripts/rebuild_atlas.py from Data/coverage_basis.csv,",
        "% Data/_rebuild/truth_dated.csv and Data/dataset_table_registry.csv",
        r"\newcommand{\DTnCanonical}{%d}" % len(canon),
        r"\newcommand{\DTnCanonicalEuk}{%d}"
        % sum(1 for r in canon if "gtdb_clusters" not in r["flags"]),
        r"\newcommand{\DTnRows}{%d}" % (len(shipped) + len(agg)),
        r"\newcommand{\DTnAggregate}{%d}" % n_agg,
        r"\newcommand{\DTnUnrepresented}{%d}" % len(unrep),
        r"\newcommand{\DTdescribedUnrepresented}{%s}"
        % thousands(sum(as_int(r["described"], 0) for r in unrep)),
        r"\newcommand{\DTdescribedEuk}{%s}" % thousands(euk_d),
    ]
    return "\n".join(lines) + "\n"


S6_ALT_HEADER = ["row_id", "partition_name", "tree_label", "citation_key",
                 "doi", "source_file", "raw_leaves", "distinct_labels",
                 "is_chronogram", "roottotip_cv", "molecular_tip_fraction",
                 "dating_code", "otl_pipeline_involvement", "archival_class",
                 "posterior_n", "imputation_status", "source_of_truth",
                 "conflict_note"]


def build_table_s6_alternative_trees(truth, warnings):
    """S6 in the manuscript body's sense: the alternative (non-canonical) trees.

    Their ultrametricity is taken from truth_dated.csv -- i.e. MEASURED from the
    shipped file -- for every alternative that is a single deposited tree."""
    rows = []
    for i, r in enumerate(sorted(truth.alt_prov,
                                 key=lambda x: (x.get("partition_name") or "",
                                                x.get("filename") or "")), 1):
        raw = r.get("partition_name") or ""
        m = re.match(r"^(.*?)\s*\(alternative:\s*(.*)\)\s*$", raw)
        part, label = (m.group(1), m.group(2)) if m else (raw, "")
        fname = (r.get("filename") or "").strip()
        d = truth.dated.get(fname)
        t = truth.tips.get(fname)
        rows.append(OrderedDict([
            ("row_id", i),
            ("partition_name", part),
            ("tree_label", label),
            ("citation_key", (r.get("citation_key") or "").strip()),
            ("doi", (r.get("doi") or "").strip()),
            ("source_file", fname),
            ("raw_leaves", as_int((t or {}).get("raw_leaves"), "")),
            ("distinct_labels", as_int((t or {}).get("distinct_labels"), "")),
            ("is_chronogram", (d or {}).get("is_chronogram", "")),
            ("roottotip_cv", fnum(as_float((d or {}).get("roottotip_cv")),
                                  FRAC_DP) if d else ""),
            ("molecular_tip_fraction",
             fnum(as_float(r.get("molecular_tip_fraction")), FRAC_DP)),
            ("dating_code", (r.get("dating_method_code") or "").strip()),
            ("otl_pipeline_involvement",
             (r.get("otl_pipeline_involvement") or "").strip()),
            ("archival_class", (r.get("archival_class") or "").strip()),
            ("posterior_n", (r.get("posterior_n") or "").strip()),
            ("imputation_status", (r.get("imputation_status") or "").strip()),
            ("source_of_truth", (r.get("source_of_truth") or "").strip()),
            ("conflict_note", (r.get("conflict_note") or "").strip()),
        ]))
    return render_csv(S6_ALT_HEADER, rows)


# -------------------------------------------------------------- Table S7 ----

S7_HEADER = ["study", "year", "estimate", "scope", "method", "is_this_study"]
S7_SECTION_RE = re.compile(
    r"\\section\{Supplementary Table S7.*?\\begin\{longtable\}(.*?)"
    r"\\end\{longtable\}", re.S)
_S7_CACHE = {}


def _s7_literature_rows(warnings):
    """The published-estimate rows of S7, carried through verbatim.

    S7 is a LITERATURE COMPILATION -- other authors' published estimates.  It is
    not derivable from the atlas's own data, so it is PASSTHROUGH; only the
    final "This study" row is recomputed.  supplementary_body.tex is opened
    READ-ONLY: another phase owns that file."""
    if "rows" not in _S7_CACHE:
        _S7_CACHE["rows"] = _parse_s7(warnings)
    return _S7_CACHE["rows"]


def _parse_s7(warnings):
    if not os.path.exists(SUPPLEMENTARY_BODY):
        warnings.append("supplement %s not found; S7's literature rows cannot "
                        "be carried through and only the computed 'This study' "
                        "row is emitted" % rel(SUPPLEMENTARY_BODY))
        return []
    with open(SUPPLEMENTARY_BODY, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = S7_SECTION_RE.search(text)
    if not m:
        warnings.append("could not locate the S7 longtable in %s; only the "
                        "computed 'This study' row is emitted"
                        % rel(SUPPLEMENTARY_BODY))
        return []
    rows = []
    for line in m.group(1).splitlines():
        line = line.strip()
        if (not line or line.startswith("\\") or "textbf{Study}" in line
                or "This study" in line or not line.endswith(r"\\")):
            continue
        # Split on UNESCAPED ampersands only: several study names contain "\&"
        # and a naive split would shred the row into six cells and drop it.
        cells = [c.strip() for c in re.split(r"(?<!\\)&", line[:-2])]
        if len(cells) != len(S7_HEADER) - 1:
            warnings.append("S7 row not in %d-column form, skipped: %s"
                            % (len(S7_HEADER) - 1, line[:70]))
            continue
        rows.append(OrderedDict(zip(S7_HEADER, cells + ["no"])))
    if not rows:
        warnings.append("S7 longtable found but no literature rows parsed from "
                        "it; check the table formatting in %s"
                        % rel(SUPPLEMENTARY_BODY))
    return rows


def _s7_this_study_row(head):
    """The one S7 row that is ours -- recomputed, never copied."""
    al = head["denominators"]["all_life"]
    million = 1e6
    lo, hi = al["estimated_low"], al["estimated_high"]
    tot = al["estimated_total"]
    return OrderedDict([
        ("study", "This study"),
        ("year", ""),      # a publication year is not a derivable quantity
        ("estimate", "%.1fM (range %.1f--%.1fM)"
         % (tot / million, lo / million, hi / million)),
        ("scope", "%d partitions (all life)" % al["n_partitions"]),
        ("method", "Per-group authoritative estimates"),
        ("is_this_study", "yes"),
    ])


def build_table_s7_csv(head, warnings):
    warnings.append(
        "S7: the 'This study' year cell is emitted BLANK. A publication year is "
        "not a quantity this script can derive and hard-coding one is "
        "forbidden; fill it from the manuscript metadata when wiring S7 in.")
    return render_csv(S7_HEADER,
                      _s7_literature_rows(warnings) + [_s7_this_study_row(head)])


def build_table_s7_tex(head, warnings):
    cols = r"p{4.2cm} c p{3cm} p{3cm} p{3.5cm}"
    header_row = (r"\textbf{Study} & \textbf{Year} & \textbf{Estimate} & "
                  r"\textbf{Scope} & \textbf{Method} \\")
    body = [" & ".join([r["study"], r["year"], r["estimate"], r["scope"],
                        r["method"]]) + r" \\"
            for r in _s7_literature_rows(warnings)]
    body.append(r"\midrule")
    t = _s7_this_study_row(head)
    body.append(" & ".join([
        r"\textbf{This study}",
        (r"\textbf{%s}" % tex_escape(t["year"])) if t["year"] else "",
        r"\textbf{$\sim$%s}" % tex_escape(t["estimate"]),
        r"\textbf{%s}" % tex_escape(t["scope"]),
        r"\textbf{%s}" % tex_escape(t["method"]),
    ]) + r" \\")
    return _longtable(cols, header_row, body,
                      caption="Published estimates of total species on Earth.")


# ==========================================================================
# HEADLINE NUMBERS -- the single source both the figures and the manuscript read
# ==========================================================================

def compute_headline(truth, basis, shipped, mol, ledger, pt, warnings):
    euk_rows = [r for r in basis if is_euk(r)]
    prok_rows = [r for r in basis if not is_euk(r)]
    den_euk = sum(as_int(r["described"], 0) for r in euk_rows)
    den_prok = sum(as_int(r["described"], 0) for r in prok_rows)
    den_all = den_euk + den_prok
    if den_euk <= 0:
        die("the eukaryotic denominator computed as %d; refusing to emit "
            "percentages against a non-positive denominator." % den_euk)

    rungs = compute_rungs(basis, shipped, truth)

    # -- per category, eukaryotic only -------------------------------------
    cats = OrderedDict()
    order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    for cat in sorted({r["category"] for r in euk_rows},
                      key=lambda c: (order.get(c, len(order)), c)):
        rows_c = [r for r in euk_rows if r["category"] == cat]
        d = sum(as_int(r["described"], 0) for r in rows_c)
        c0 = sum(as_int(r["c0_units"], 0) for r in rows_c)
        c1 = sum(as_int(r["c0_units"], 0) * rung_weight(r, "C1")
                 for r in rows_c)
        c1d = sum(as_int(r["c0_units"], 0) * rung_weight(r, "C1d")
                  for r in rows_c)
        # REACHABILITY.  C1 and C1d are weighted sums and are genuinely
        # fractional; publishing only the ROUNDED integer beside an unrounded
        # percentage leaves a percentage that no printed pair reproduces --
        # which is exactly the defect 23.7 % was.  Both the exact sum and the
        # rounded one are therefore published, as the top-level rungs do.
        cats[cat] = OrderedDict([
            ("described", d),
            ("denominator", d),
            ("n_partitions", len(rows_c)),
            ("C0_sum", c0), ("C0_pct", pct(c0, d)),
            ("C1_sum_exact", round(c1, NUMERATOR_DP)),
            ("C1_sum", int(round(c1))), ("C1_pct", pct(c1, d)),
            ("C1d_sum_exact", round(c1d, NUMERATOR_DP)),
            ("C1d_sum", int(round(c1d))), ("C1d_pct", pct(c1d, d)),
        ])

    # -- prokaryotes, on their own line, never mixed ------------------------
    prok = OrderedDict([
        ("n_partitions", len(prok_rows)),
        ("described", den_prok),
        ("note", "reported on its own line. Prokaryotes NEVER enter a "
                 "eukaryotic numerator or denominator."),
    ])
    p_named = sum(as_int(r["c0_units"], 0) for r in prok_rows)
    prok["C0_sum"] = p_named
    prok["C0_pct"] = pct(p_named, den_prok)
    prok["per_partition"] = OrderedDict()
    for r in sorted(prok_rows, key=lambda x: x["partition"]):
        d = as_int(r["described"], 0)
        prok["per_partition"][r["partition"]] = OrderedDict([
            ("described", d),
            ("c0_units", as_int(r["c0_units"], 0)),
            ("C0_pct", pct(as_int(r["c0_units"], 0), d)),
            ("is_chronogram", r["is_chronogram"]),
            ("roottotip_cv", as_float(r["roottotip_cv"])),
            ("C1d_sum", int(round(as_int(r["c0_units"], 0)
                                  * rung_weight(r, "C1d")))),
        ])

    # -- tree structure ------------------------------------------------------
    n_trees = len(truth.dated)
    n_chrono = sum(1 for r in truth.dated.values() if truthy(r["is_chronogram"]))
    canon_chrono = [p for p in truth.canon_prov if truth.is_chronogram(p)]
    trees = OrderedDict([
        ("n_shipped", n_trees),
        ("n_chronograms_all_shipped", n_chrono),
        ("n_partition_canonical", len(truth.canon_prov)),
        ("n_partition_canonical_eukaryotic", len(truth.canon_euk)),
        ("n_partition_canonical_chronograms", len(canon_chrono)),
        ("chronogram_partitions", sorted(canon_chrono)),
        ("relative_time_partitions", sorted(
            p for p in truth.canon_prov
            if truth.is_geometric_chronogram(p)
            and truth.time_scale(p) == "relative")),
        ("rule", "root-to-tip coefficient of variation measured on the SHIPPED "
                 ".nwk file with the audit's recursive parser (truth_dated.csv), "
                 "AND branch lengths in absolute time (Data/time_scale_checks.csv; "
                 "D14: an ultrametric tree in relative time is its own class); "
                 "the ledger's `dated` column is never read"),
        ("note", "everything at C1d and below rests on these trees ONLY"),
    ])

    archival = defaultdict(int)
    for part in truth.canon_prov:
        archival[truth.archival_class(part)] += 1

    disc = disclosures(basis, shipped, truth, mol, rungs, warnings)

    head = OrderedDict([
        ("_README", (
            "Generated by scripts/rebuild_atlas.py from Data/_rebuild/truth_*.csv "
            "and the shipped trees. EVERY value here is computed at run time; "
            "none is hard-coded. This file is the single source for the "
            "manuscript's and the figures' headline quantities -- read it, do "
            "not retype the numbers. Numerator and denominator are column sums "
            "of Data/coverage_basis.csv, so a denominator that exists in no "
            "file cannot be written down.")),
        ("generator", "scripts/rebuild_atlas.py"),
        ("basis_file", rel(OUT["coverage_basis.csv"])),
        ("truth_files", [rel(p) for p in TRUTH_FILES.values()]),
        ("rule_parameters", OrderedDict([
            ("source_numerator_tolerance", SOURCE_NUMERATOR_TOL),
            ("source_nonsubset_tolerance", SOURCE_NONSUBSET_TOL),
            ("source_scan_max_bytes", SOURCE_SCAN_MAX_BYTES),
            ("cv_recheck_tolerance", CV_RECHECK_TOL),
            ("chronogram_tolerance", CHRONOGRAM_TOLERANCE),
            ("elasticity_fraction", ELASTICITY_FRACTION),
            ("note", "RULES, not reported quantities. No coverage figure "
                     "anywhere in this file is a literal."),
        ])),
        ("scope", OrderedDict([
            ("numerator_universe",
             "the %d EUKARYOTIC CANONICAL trees" % len(truth.canon_euk)),
            ("denominator_universe",
             "described species over the %d eukaryotic partitions"
             % len(euk_rows)),
            ("prokaryotes", "reported separately; never mixed"),
            ("all_life", "SI only, and Figure 1, always labelled as such"),
        ])),
        ("denominators", OrderedDict([
            ("eukaryotic", OrderedDict([
                ("n_partitions", len(euk_rows)),
                ("described", den_euk),
                ("note", "THE denominator for every rung."),
            ])),
            ("prokaryotic", OrderedDict([
                ("n_partitions", len(prok_rows)),
                ("described", den_prok),
            ])),
            ("all_life", OrderedDict([
                ("n_partitions", len(basis)),
                ("described", den_all),
                ("estimated_total", sum(pt.estimated_total(r["partition"]) or 0
                                        for r in basis)),
                ("estimated_low",
                 sum(as_int(truth.species[r["partition"]].get("estimated_low"), 0)
                     for r in basis)),
                ("estimated_high",
                 sum(as_int(truth.species[r["partition"]].get("estimated_high"), 0)
                     for r in basis)),
                ("note", "SI and Figure 1 ONLY, and it must be labelled "
                         "all-life there. NEVER a coverage denominator."),
            ])),
        ])),
        ("rungs", rungs),
        ("rung_definitions", OrderedDict(RUNGS)),
        ("naming_discipline", (
            "C1 is MOLECULARLY SUPPORTED PLACEMENT, operationalised by tip "
            "evidence within the analysis that produced the tree. It is NOT a "
            "predicate on placement traceability: f_mol is an axis-(a) "
            "quantity. The genuine axis-(b) distinction is reported as the "
            "NOMINAL annotation topology_assembly and never as a rung, because "
            "the deposit carries no measured per-tip placement-warrant flag. "
            "Methods must say so plainly; claiming the two-axis conflict was "
            "'dissolved' is an overclaim a referee will catch.")),
        ("monotonicity", (
            "the rungs are nested weighted sums with weights in [0,1], hence "
            "non-increasing; verified per partition, %d partitions, zero "
            "violations. NESTED TIP SETS ARE NOT CLAIMED: the computed object "
            "is an expectation, and the advertised proof describes the "
            "computed object." % len(euk_rows))),
        ("units", (
            "C0 counts SPECIES. C1 and below count TIP-EQUIVALENTS wherever the "
            "partition is not enumerated in %s. Do not write 'species' at C1 "
            "or below for a partition that is not listed there."
            % rel(OUT["molecular_tips.csv"]))),
        ("per_category_eukaryotic", cats),
        ("prokaryotic_line", prok),
        ("trees", trees),
        ("archival_classes_canonical",
         OrderedDict((k, archival[k]) for k in sorted(archival) if k)),
        ("topology_assembly", OrderedDict([
            ("classifier", "automated_synthesis if otl_pipeline_involvement == "
                           "'%s'; compiled if inference in %s; else "
                           "primary_analysis"
             % (ASSEMBLY_OTL_TRIGGER, list(ASSEMBLY_INFERENCE_TRIGGERS))),
            ("source", "derived from UNRECONCILED fields of %s (the pre-audit "
                       "ledger). Methods MUST state this. It is an ANNOTATION, "
                       "not a criterion: nominal, unordered, no strong/weak "
                       "reading." % rel(TREE_PROVENANCE)),
            ("counts", OrderedDict(
                (c, sum(1 for r in euk_rows if r["topology_assembly"] == c))
                for c in sorted({r["topology_assembly"] for r in euk_rows
                                 if r["topology_assembly"]}))),
        ])),
        ("sensitivities", OrderedDict([
            ("f_mol_bands", band_sensitivity(basis, shipped, truth, warnings)),
            ("leave_assembly_class_out",
             assembly_sensitivity(basis, shipped, truth)),
            ("denominator", denominator_sensitivity(basis, shipped, truth)),
            ("c1_definition", c1_definition_sensitivity(basis, shipped, truth)),
            ("c1b_calibration_fixed",
             c1b_calibration_sensitivity(basis, shipped, truth)),
            ("c1b_all_intervals",
             c1b_all_intervals_sensitivity(basis, shipped, truth)),
            ("c1b_crown_calibrated",
             c1b_crown_calibrated_sensitivity(basis, shipped, truth)),
            ("c1b_rule_p_counted",
             c1b_rule_p_sensitivity(basis, shipped, truth)),
            ("c0_gbif_unmatched", c0_gbif_sensitivity(basis, shipped)),
        ])),
        ("disclosures", disc),
        ("molecular_evidence", OrderedDict(
            (p, OrderedDict([
                ("f_mol", mol[p].f_mol),
                ("f_basis", mol[p].f_basis),
                ("f_low", mol[p].f_low),
                ("f_high", mol[p].f_high),
                ("molecular_source_file", mol[p].source_file),
                ("measured_intersection", mol[p].measured_intersection),
                ("exact_label_intersection", mol[p].exact_label_intersection),
                ("recorded_numerator", mol[p].recorded_numerator),
                ("contested_source", mol[p].contested_source),
                ("contested_f_mol", mol[p].contested_f),
                ("notes", mol[p].notes),
            ])) for p in sorted(mol))),
        ("tip_convention",
         "distinct standardised labels in the shipped Newick file "
         "(truth_tips.csv `tips_adopted`)"),
        ("retired", RETIRED),
        ("open_flags", OPEN_FLAGS),
        ("warnings", warnings),
    ])
    return head


def name_resolution(truth, shipped, basis):
    """S05: per eukaryotic canonical partition, how many distinct shipped
    labels resolve to a Catalogue of Life accepted name (Data/standardized/
    col_resolution.csv, written by scripts/match_col.py against the COL26.9
    base release, D18), how many do not, and how many distinct accepted names
    they reach.  The distinct accepted names ARE the partition's C0 count
    (basis column c0_units); prokaryotic partitions are counted in GTDB species
    clusters instead (D12) and carry no CoL row."""
    euk = {r["partition"] for r in basis if is_euk(r)}
    den_euk = sum(as_int(r["described"], 0) for r in basis if is_euk(r))
    per = OrderedDict()
    tot_labels = tot_resolved = 0
    accepted_union = set()
    tot_accepted = 0
    for part in sorted(truth.canon_prov):
        if part not in euk:
            continue
        rows = shipped.col_rows.get(part, {})
        types = defaultdict(int)
        accepted = set()
        resolved = 0
        for sid in sorted(set(shipped.ids[part].values())):
            r = rows[sid]
            types[r["match_type"]] += 1
            if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
                resolved += 1
                accepted.add(r["col_id"])
        n = len(shipped.ids[part])
        per[part] = OrderedDict([
            ("eukaryotic", True),
            ("labels", n),
            ("resolved_labels", resolved),
            ("distinct_accepted_names", len(accepted)),
            ("unresolved_share", OrderedDict([
                ("numerator", n - resolved), ("denominator", n),
                ("percent", pct(n - resolved, n))])),
            ("match_types", OrderedDict(sorted(types.items()))),
        ])
        tot_labels += n
        tot_resolved += resolved
        tot_accepted += len(accepted)
        accepted_union |= accepted
    return OrderedDict([
        ("authority", "Catalogue of Life COL26.9 base release (ChecklistBank "
                      "316321, doi:10.48580/dgz5n), through "
                      "Data/standardized/col_resolution.csv"),
        ("resolved_match_types", list(COL_RESOLUTION_RESOLVED)),
        ("rule", "a tip counts at C0 when it resolves to a CoL accepted name; "
                 "tips resolving to the same accepted name count once"),
        ("eukaryotic", OrderedDict([
            ("unresolved_share", OrderedDict([
                ("numerator", tot_labels - tot_resolved),
                ("denominator", tot_labels),
                ("percent", pct(tot_labels - tot_resolved, tot_labels))])),
            ("C0_resolved_labels", OrderedDict([
                ("numerator", tot_resolved), ("denominator", den_euk),
                ("percent", pct(tot_resolved, den_euk))])),
            ("C0_distinct_accepted_names_sum", OrderedDict([
                ("numerator", tot_accepted), ("denominator", den_euk),
                ("percent", pct(tot_accepted, den_euk))])),
            ("C0_distinct_accepted_names_union", OrderedDict([
                ("numerator", len(accepted_union)), ("denominator", den_euk),
                ("percent", pct(len(accepted_union), den_euk))])),
        ])),
        ("per_partition", per),
    ])


def build_c0_name_resolution(head, warnings):
    rows = []
    for part, v in head["name_resolution_S05"]["per_partition"].items():
        rows.append(OrderedDict([
            ("partition", part),
            ("eukaryotic", "yes" if v["eukaryotic"] else "no"),
            ("labels", v["labels"]),
            ("resolved_labels", v["resolved_labels"]),
            ("unresolved_labels", v["unresolved_share"]["numerator"]),
            ("unresolved_pct", fnum(v["unresolved_share"]["percent"], PCT_DP)),
            ("distinct_accepted_names", v["distinct_accepted_names"]),
            ("match_types", "; ".join("%s=%d" % kv
                                      for kv in v["match_types"].items())),
        ]))
    return render_csv(["partition", "eukaryotic", "labels", "resolved_labels",
                       "unresolved_labels", "unresolved_pct",
                       "distinct_accepted_names", "match_types"], rows)


def revision_checks(truth, shipped, names, warnings):
    """D7 (off-target tips) and D8 (deposited chronograms), read from their
    registers and measured against the shipped files.

    D7: Data/offtarget_checks.csv records which canonical trees have been
    checked for tips outside their partition; Data/pruned_tips.csv lists every
    tip pruned from a canonical tree.  A pruned id still present in its shipped
    tree is a refusal.  The share of canonical trees not yet checked is
    published, so an unchecked tree is visible as unchecked.

    D8: Data/dating_deposit_checks.csv names, per partition, the chronogram the
    source study deposited.  How many shipped tips it covers is MEASURED when
    the file is present; the outcome follows from that count, never from a
    typed verdict."""
    checked = {(r.get("filename") or "").strip(): r
               for r in read_csv(OFFTARGET_CHECKS, required=False)}
    pruned = defaultdict(list)
    for r in read_csv(PRUNED_TIPS, required=False):
        pruned[(r.get("filename") or "").strip()].append(
            (r.get("standardized_id") or "").strip())
    per = OrderedDict()
    for part in sorted(truth.canon_prov):
        fn = truth.canon_file[part]
        ids = pruned.get(fn, [])
        still = sorted(i for i in ids
                       if names.by_id.get(i, i) in shipped.labels[part])
        if still:
            die("%s: %d tip(s) listed in %s are still in the shipped tree "
                "(e.g. %s)" % (fn, len(still), rel(PRUNED_TIPS), still[:5]))
        if ids and fn not in checked:
            die("%s has pruned tips but no row in %s"
                % (fn, rel(OFFTARGET_CHECKS)))
        per[part] = OrderedDict([("canonical_tree", fn),
                                 ("checked", fn in checked),
                                 ("n_pruned", len(ids))])
    unchecked = [p for p, v in per.items() if not v["checked"]]
    offtarget = OrderedDict([
        ("rule", "D7: tips outside the partition are pruned from the canonical "
                 "tree before counting and before crown-node extraction"),
        ("n_canonical", len(per)),
        ("n_checked", len(per) - len(unchecked)),
        ("n_pruned_tips", sum(v["n_pruned"] for v in per.values())),
        ("unchecked_partitions", unchecked),
        ("per_partition", per),
    ])

    dating = OrderedDict()
    for r in read_csv(DATING_DEPOSIT_CHECKS, required=False):
        part = (r.get("partition") or "").strip()
        if part not in truth.canon_prov:
            die("%s names %r, which is not a canonical partition"
                % (rel(DATING_DEPOSIT_CHECKS), part))
        rel_src = (r.get("deposited_chronogram") or "").strip()
        src = os.path.join(ROOT, "Data", rel_src)
        n = truth.labels_unique(part)
        covered = None
        if not rel_src:
            outcome = "none_deposited"
        elif os.path.exists(src):
            covered = len(_source_sets(src, names, truth.group(part))[0]
                          & shipped.labels[part])
            outcome = ("chronogram_covers_shipped_tips" if covered == n
                       else "gap_in_deposit")
        else:
            outcome = "not_rechecked_source_absent"
            warnings.append("%s: %s is not present; D8 coverage not measured "
                            "on this run" % (part, rel_src))
        dating[part] = OrderedDict([
            ("canonical_tree", truth.canon_file[part]),
            ("deposited_chronogram", rel_src),
            ("shipped_tips", n),
            ("chronogram_tips_in_shipped", covered),
            ("outcome", outcome),
        ])
    return OrderedDict([("offtarget_D7", offtarget),
                        ("deposited_chronogram_D8", OrderedDict([
                            ("rule", "D8: where the same study deposited a "
                                     "chronogram covering at least the shipped "
                                     "tips, it becomes the canonical file; "
                                     "otherwise the gap is reported"),
                            ("per_partition", dating)]))])


CROWN_HEADER = ["partition", "canonical_tree", "crown_clade", "crown_age_ma",
                "age_statistic", "geometry_root_ma", "interval_kind",
                "counts_for_c1b", "c1b_exclusion", "interval_lo_ma",
                "interval_center_ma", "interval_hi_ma", "center_statistic",
                "n_trees", "interval_width_myr", "degenerate", "at_bound",
                "calibration", "prior_reproduced", "interval_source"]
CROWN_CLADES = os.path.join(ROOT, "Data", "crown_clades.csv")
# kinds whose centre is the deposit's own summary of the node age (D20)
CROWN_CENTER_KINDS = ("node_hpd_annotation", "posterior_sample",
                      "bootstrap_dating_sample", "calibration_prior_sample")
CROWN_SUMMARY_KINDS = ("node_hpd_annotation", "node_label_interval")


def crown_age_rule(iv, root_h):
    """(crown age, what it is) by decision D20.  The crown age is the deposit's
    own summary of the focal node: the annotated median height of a
    TreeAnnotator summary tree, or the median height over a dated tree set
    (TreeAnnotator's median-heights convention); otherwise the focal node's
    height in the deposited chronogram; with no interval register row, the
    root height of the canonical file.  The branch-length geometry of a
    summary tree is kept as geometry_root_ma, a diagnostic."""
    if iv and iv["kind"] in CROWN_CENTER_KINDS and iv["interval_center_ma"]:
        stat = ("annotated %s" % iv["center_statistic"]
                if iv["kind"] == "node_hpd_annotation"
                else "median over %s trees" % iv["n_trees"])
        return as_float(iv["interval_center_ma"]), stat
    if iv and iv["kind"] in CROWN_SUMMARY_KINDS:
        return (as_float(iv["source_point_height_ma"]),
                "focal node height in the deposited chronogram")
    return root_h, "root height of the canonical file"


def crown_ages(truth, basis, warnings):
    """S02.  One row per eukaryotic partition whose canonical file is a
    chronogram in absolute time.

    The focal clade is named in Data/crown_clades.csv: after D7 it is the clade
    of the retained tips, or, where the source's ingroup is a subclade whose
    outgroups belong to the partition, that ingroup (D21).  crown_age_ma
    follows crown_age_rule() (D20).  For a summary-tree deposit the focal node
    is identified by geometry: with no kept outgroup its height in the deposit
    equals the canonical root height; with kept outgroups the deposit's root
    equals the canonical root and the focal node lies below it."""
    reg = crown_intervals()
    clades = {r["partition"]: r for r in read_csv(CROWN_CLADES)}
    rows = []
    for r in basis:
        part = r["partition"]
        if not is_euk(r) or r["time_scale"] != "absolute" \
                or r["is_chronogram"] != "yes":
            continue
        fn = r["canonical_tree"]
        h = ROOT_HEIGHT.get(fn)
        if h is None:
            die("no root height measured for %s" % fn)
        if not (clades.get(part) or {}).get("crown_clade"):
            die("%s: no crown_clade in %s" % (part, rel(CROWN_CLADES)))
        iv = reg.get(part)
        age, stat = crown_age_rule(iv, h)
        row = OrderedDict((k, "") for k in CROWN_HEADER)
        row.update(partition=part, canonical_tree=fn,
                   crown_clade=clades[part]["crown_clade"],
                   crown_age_ma=fnum(age, CROWN_DP), age_statistic=stat,
                   geometry_root_ma=fnum(h, CROWN_DP))
        if iv:
            if iv["kind"] in CROWN_SUMMARY_KINDS:
                sp = as_float(iv["source_point_height_ma"])
                if as_int(iv.get("n_focal_outgroups"), 0) == 0:
                    if sp is None or abs(sp - h) > CROWN_POINT_TOL:
                        die("%s: the focal node of %s sits at %s Ma but the "
                            "canonical root is at %.6f Ma; the register's "
                            "'outside' list does not reproduce the D7 pruning"
                            % (part, iv["source"],
                               iv["source_point_height_ma"], h))
                else:
                    sr = as_float(iv["source_root_height_ma"])
                    if sr is None or abs(sr - h) > CROWN_POINT_TOL or sp >= h:
                        die("%s: with kept outgroups the deposit's root (%s) "
                            "must equal the canonical root (%.6f) and the "
                            "focal node (%s) lie below it"
                            % (part, sr, h, sp))
            row.update(
                interval_kind=iv["kind"], counts_for_c1b=iv["counts_for_c1b"],
                c1b_exclusion=iv["c1b_exclusion"],
                interval_lo_ma=fnum(as_float(iv["interval_lo_ma"]), CROWN_DP),
                interval_center_ma=(fnum(as_float(iv["interval_center_ma"]),
                                         CROWN_DP)
                                    if iv["interval_center_ma"] else ""),
                interval_hi_ma=fnum(as_float(iv["interval_hi_ma"]), CROWN_DP),
                center_statistic=iv["center_statistic"],
                n_trees=iv["n_trees"],
                interval_width_myr=fnum(as_float(iv["interval_width_myr"]),
                                        CROWN_DP),
                degenerate=iv["degenerate"], at_bound=iv["at_bound"],
                calibration=iv.get("calibration", ""),
                prior_reproduced=iv.get("prior_reproduced", ""),
                interval_source=iv["source"])
        rows.append(row)
    rows.sort(key=lambda x: x["partition"])
    for part in reg:
        if part not in {x["partition"] for x in rows}:
            warnings.append("%s has a row in %s but its canonical tree is not "
                            "an absolute chronogram; no crown age is reported"
                            % (part, rel(CROWN_INTERVALS)))
    counted = [x for x in rows if x["counts_for_c1b"] == "yes"]
    excl = Counter(x["c1b_exclusion"] for x in rows if x["c1b_exclusion"])
    return OrderedDict([
        ("definition",
         "crown age of the focal clade (Data/crown_clades.csv, D21) by "
         "decision D20: the deposit's annotated median height, or the median "
         "height over its dated tree set, else the focal node's height in the "
         "deposited chronogram, else the canonical root height; the interval "
         "is the one the deposit carries (Data/crown_intervals.csv, "
         "scripts/crown_intervals.py). An interval counts for C1b when its "
         "kind counts and it is neither degenerate nor at a calibration bound "
         "(D22) nor a reproduction of the focal node's calibration prior "
         "(D25, Data/calibration_register.csv)."),
        ("n_absolute_chronograms", len(rows)),
        ("n_with_counted_interval", len(counted)),
        ("n_degenerate_counted_intervals",
         sum(1 for x in counted if x["degenerate"] == "yes")),
        ("n_excluded_from_c1b", OrderedDict(sorted(excl.items()))),
        ("n_age_differs_from_geometry",
         sum(1 for x in rows if abs(as_float(x["crown_age_ma"])
                                    - as_float(x["geometry_root_ma"]))
             > CROWN_POINT_TOL)),
        ("n_crown_inside_interval",
         sum(1 for x in counted
             if as_float(x["interval_lo_ma"]) <= as_float(x["crown_age_ma"])
             <= as_float(x["interval_hi_ma"]))),
        ("per_partition", OrderedDict((x["partition"], OrderedDict(
            (k, x[k]) for k in CROWN_HEADER if k != "partition"))
            for x in rows)),
    ])


def build_crown_ages_csv(head, warnings):
    rows = [OrderedDict([("partition", p)] + list(v.items()))
            for p, v in head["crown_ages"]["per_partition"].items()]
    return render_csv(CROWN_HEADER, rows)


# interval kinds that come from a set of trees (the others are intervals on
# one summary tree)
CROWN_SET_KINDS = ("posterior_sample", "bootstrap_dating_sample",
                   "calibration_prior_sample", "taxon_addition_sample",
                   "date_assignment_sample")


def crown_macro_name(partition):
    """'Seed plants' -> 'SeedPlants': letters only, so it is a TeX name."""
    return "".join(w[:1].upper() + w[1:] for w in re.split(r"[^A-Za-z]+",
                                                         partition) if w)


def lr_mixture_p(lr):
    """P(T >= lr) for T the 50:50 mixture of a point mass at zero and
    chi2(1): 0.5 * P(chi2(1) >= lr) for lr > 0."""
    if lr <= 0:
        return 1.0
    return 0.5 * math.erfc(math.sqrt(lr / 2.0))


def p_fmt(p):
    """A p-value for math mode: three decimals, or a power of ten below 0.001."""
    if p >= 0.001:
        return "%.3f" % p
    e = int(math.floor(math.log10(p)))
    return "%.1f \\times 10^{%d}" % (p / 10 ** e, e)


def prior_q_tex(q):
    """A prior quantile in prose: two decimals, or '$>$0.975' beyond a soft
    bound."""
    if q[:1] in "<>":
        return "$%s$%s" % (q[0], q[1:])
    return "%.2f" % float(q)


def build_crown_ages_macros(head, warnings):
    """Audit A13.  Every crown age, interval and width the prose or a caption
    quotes is a macro read from the one register (crown_ages.csv), so no
    second register can drift from it.  Values to one decimal; the width as a
    percentage of the crown age."""
    out = ["% generated by scripts/rebuild_atlas.py from "
           "Paper/inputs/crown_ages.csv (decisions D20-D22, D25)"]
    for part, v in head["crown_ages"]["per_partition"].items():
        n = crown_macro_name(part)
        age = as_float(v["crown_age_ma"])
        out.append("\\newcommand{\\CrownAge%s}{%.1f}" % (n, age))
        lo, hi = as_float(v["interval_lo_ma"]), as_float(v["interval_hi_ma"])
        if lo is not None and hi is not None:
            out.append("\\newcommand{\\CrownLo%s}{%.1f}" % (n, lo))
            out.append("\\newcommand{\\CrownHi%s}{%.1f}" % (n, hi))
            out.append("\\newcommand{\\CrownWidth%s}{%.1f}" % (n, hi - lo))
            out.append("\\newcommand{\\CrownWidthPct%s}{%.1f}"
                       % (n, 100.0 * (hi - lo) / age))
    # round F3 (F3-14, F3-11): where each interval with a stated prior
    # density lies in it, and the D25 measures of every two-sided prior
    for part, iv in sorted(crown_intervals().items()):
        n = crown_macro_name(part)
        if iv.get("prior_q_center"):
            for k, nm in (("prior_q_lo", "Lo"), ("prior_q_center", "Center"),
                          ("prior_q_hi", "Hi")):
                out.append("\\newcommand{\\CrownPriorQ%s%s}{%s}"
                           % (nm, n, prior_q_tex(iv[k])))
        if iv.get("prior_containment"):
            out.append("\\newcommand{\\CrownPriorIn%s}{%.0f}"
                       % (n, 100 * as_float(iv["prior_containment"])))
            out.append("\\newcommand{\\CrownPriorRatio%s}{%.2f}"
                       % (n, as_float(iv["prior_width_ratio"])))
            out.append("\\newcommand{\\CrownPriorLo%s}{%.1f}"
                       % (n, as_float(iv["prior_lo_ma"])))
            out.append("\\newcommand{\\CrownPriorHi%s}{%.1f}"
                       % (n, as_float(iv["prior_hi_ma"])))
    # round F2: the range of widths of the usable intervals
    wc = [100.0 * (as_float(v["interval_hi_ma"]) - as_float(v["interval_lo_ma"]))
          / as_float(v["crown_age_ma"])
          for v in head["crown_ages"]["per_partition"].values()
          if v["counts_for_c1b"] == "yes"]
    out.append("\\newcommand{\\CrownWidthPctMinCounted}{%.1f}" % min(wc))
    out.append("\\newcommand{\\CrownWidthPctMaxCounted}{%.1f}" % max(wc))
    return "\n".join(out) + "\n"


# ==========================================================================
# NUMBER MACROS -- every quantity the manuscript prose quotes
# ==========================================================================

def macro_pct(x):
    """Percent formatting rule for the prose (D9): one decimal; values below 1
    get two decimals so that a sub-percent figure does not round to 0.x."""
    if x is None:
        return "--"
    return "%.1f" % x if abs(x) >= 1.0 else "%.2f" % x


def macro_int(n):
    return thousands(int(round(float(n))))


def build_number_macros(head, basis, warnings):
    """Paper/inputs/numbers_macros.tex.  The manuscript writes no number by
    hand (AGENTS.md rule 1): every count, percentage and estimate the prose
    quotes is a macro read from headline_numbers.json or the basis table.
    Names are letters only; a partition's name is crown_macro_name()."""
    m = OrderedDict()
    den = head["denominators"]
    tr = head["trees"]
    m["NumPartitionsAll"] = macro_int(den["all_life"]["n_partitions"])
    m["NumPartitionsEuk"] = macro_int(den["eukaryotic"]["n_partitions"])
    m["NumPartitionsProk"] = macro_int(den["prokaryotic"]["n_partitions"])
    m["NumTreesShipped"] = macro_int(tr["n_shipped"])
    m["NumTreesChron"] = macro_int(tr["n_chronograms_all_shipped"])
    m["NumCanonical"] = macro_int(tr["n_partition_canonical"])
    m["NumCanonicalEuk"] = macro_int(tr["n_partition_canonical_eukaryotic"])
    m["NumCanonicalChron"] = macro_int(tr["n_partition_canonical_chronograms"])
    m["NumCanonicalRelTime"] = macro_int(len(tr["relative_time_partitions"]))
    m["DescribedEuk"] = macro_int(den["eukaryotic"]["described"])
    m["DescribedProk"] = macro_int(den["prokaryotic"]["described"])
    rep = head["sensitivities"]["denominator"]["represented_partitions_only"]
    m["DescribedRepresented"] = macro_int(rep["described"])

    # the four rungs, on the headline (union) convention
    rg = head["rungs"]
    names = (("C0", "CZero"), ("C1", "COne"), ("C1d", "COneD"),
             ("C1b", "COneB"))
    for key, nm in names:
        m[nm + "N"] = macro_int(rg[key]["numerator_union"])
        m[nm + "Pct"] = macro_pct(rg[key]["percent_union"])
        m[nm + "NPart"] = macro_int(rg[key]["n_partitions_contributing"])
    m["COneNEnumerated"] = macro_int(rg["C1"]["c1_bases"]["enumerated"])
    m["COneNAllTips"] = macro_int(rg["C1"]["c1_bases"]["all_tips"])
    asserted = sum(as_int(r["c0_units"], 0) * rung_weight(r, "C1")
                   for r in basis if is_euk(r) and r.get("c1_basis") == "all_tips")
    m["COneAssertedN"] = macro_int(asserted)
    m["COneAssertedPct"] = macro_pct(pct(asserted, rg["C1"]["numerator_sum"]))
    unrep = sum(as_int(r["described"], 0) for r in basis
                if is_euk(r) and as_int(r.get("c0_units"), 0) <= 0)
    m["UnrepresentedPct"] = macro_pct(pct(unrep, den["eukaryotic"]["described"]))
    for (a, an), (b, bn) in zip(names[1:], names[:-1]):
        m[an + "Of" + bn] = macro_pct(pct(rg[a]["numerator_union"],
                                          rg[b]["numerator_union"]))

    dom = head["disclosures"]["dominance"]
    for key, nm in names:
        m["SeedPlantsShare" + nm] = macro_pct(dom["share_of_%s_pct" % key])
    # the eukaryote-wide levels with the dominant group's tree and its
    # described species both removed (Limitations)
    dom_d = sum(as_int(r["described"], 0) for r in basis
                if r["partition"] == dom["partition"])
    for key, nm in names[:3]:
        m["Excl" + crown_macro_name(dom["partition"]) + nm + "Pct"] = macro_pct(
            pct(dom["%s_total" % key] - dom["contribution_to_%s" % key],
                den["eukaryotic"]["described"] - dom_d))

    for cat, v in head["per_category_eukaryotic"].items():
        cn = crown_macro_name(cat)
        for key, nm in names[:3]:
            m["Cat%s%sPct" % (cn, nm)] = macro_pct(v["%s_pct" % key])

    # R1-11: tips that do not count as a named species (unresolved labels and
    # tips merged into one accepted name), and named coverage counted in tips
    canon_e = [r for r in basis if is_euk(r) and as_int(r.get("c0_units"), 0) > 0]
    m["TipsNotNamedN"] = macro_int(sum(as_int(r["labels_unique"], 0)
                                       - as_int(r["c0_units"], 0) for r in canon_e))
    m["SensTipsCZeroPct"] = macro_pct(pct(sum(
        min(as_int(r["labels_unique"], 0), as_int(r["described"], 0))
        for r in canon_e), den["eukaryotic"]["described"]))
    for r in canon_e:
        m["TipsNotNamed" + crown_macro_name(r["partition"])] = macro_int(
            as_int(r["labels_unique"], 0) - as_int(r["c0_units"], 0))
    # R1-06: groups whose reference tree holds taxonomically placed tips, and
    # how many of their sources distribute a tree of the molecular species
    mol_n = _molecular_tip_counts()
    taxo = sorted(r["partition"] for r in canon_e
                  if r["c1_basis"] == "enumerated"
                  and mol_n[r["partition"]] < as_int(r["labels_unique"], 0))
    molonly = _molecular_only_trees()
    if sorted(molonly) != taxo:
        die("%s lists %s; the groups with taxonomically placed tips are %s"
            % (rel(MOLECULAR_ONLY_TREES), sorted(molonly), taxo))
    m["NGroupsTaxonomicTips"] = macro_int(len(taxo))
    m["NMolOnlyTrees"] = macro_int(len(molonly))
    # F3-16: which of those trees are dated, undated phylograms, or carry no
    # branch lengths, measured on the distributed files
    mod = _molecular_only_dating()
    m["NMolOnlyDated"] = macro_int(sum(v == "time" for v in mod.values()))
    m["MolOnlyDatedNames"] = name_list(sorted(p for p, v in mod.items()
                                              if v == "time"))
    m["MolOnlyPhylogramNames"] = name_list(sorted(
        p for p, v in mod.items() if v == "substitutions"))
    m["MolOnlyNoLengthNames"] = name_list(sorted(
        p for p, v in mod.items() if v == "none"))
    # R1-09: what the study behind each reference tree that is not an
    # absolute chronogram did about dating
    sd = _study_dating()
    for key, nm in STUDY_DATING_MACROS.items():
        m["StudyDating%sN" % nm] = macro_int(sum(
            1 for r in canon_e if r["partition"] in sd and sd[r["partition"]] == key
            and not (r["is_chronogram"] == "yes" and r["time_scale"] == "absolute")))

    nr = head["name_resolution_S05"]["eukaryotic"]["unresolved_share"]
    m["LabelsEukN"] = macro_int(nr["denominator"])
    m["UnresolvedLabelsN"] = macro_int(nr["numerator"])
    m["UnresolvedLabelsPct"] = macro_pct(nr["percent"])

    s = head["sensitivities"]
    g = s["c0_gbif_unmatched"]
    m["SensGbifSpeciesAdded"] = macro_int(g["species_added"])
    m["SensGbifCZeroPct"] = macro_pct(g["percent"])
    col = s["denominator"]["col_single_authority"]["filled"]
    m["SensColDenominator"] = macro_int(col["denominator"])
    for key, nm in names:
        m["SensCol%sPct" % nm] = macro_pct(col["rungs"][key]["percent_union"]
                                           if "percent_union" in col["rungs"][key]
                                           else col["rungs"][key]["percent"])
        m["SensRep%sPct" % nm] = macro_pct(rep["rungs"][key]["percent_union"])
    u = head["union_bound_S09"]
    m["SensUnionNTrees"] = macro_int(u["n_noncanonical_trees"])
    m["SensUnionNamesAdded"] = macro_int(u["names_added"])
    m["SensUnionCZeroPct"] = macro_pct(u["C0_union_pct"])
    m["SensUnionCOneHighPct"] = macro_pct(u["C1_union_high_pct"])
    prod = s["c1_definition"]["product"]
    m["SensProductCOnePct"] = macro_pct(prod["C1"]["percent"])
    m["SensProductCOneDPct"] = macro_pct(prod["C1d"]["percent"])
    cal = s["c1b_calibration_fixed"]["rungs"]["C1b"]
    alli = s["c1b_all_intervals"]["rungs"]["C1b"]
    m["SensAllIntCOneBN"] = macro_int(alli["numerator_union"])
    m["SensAllIntCOneBPct"] = macro_pct(alli["percent_union"])
    m["SensCalibCOneBN"] = macro_int(cal["numerator_union"])
    m["SensCalibCOneBPct"] = macro_pct(cal["percent_union"])
    cc = s["c1b_crown_calibrated"]
    bycal = cc["counted_by_crown_calibration"]
    if bycal["not_checked"]:
        die("a counted interval sits on an unchecked crown calibration "
            "(%s); the text states that none does"
            % ", ".join(bycal["not_checked"]))
    m["SensCrownCalCOneBN"] = macro_int(cc["rungs"]["C1b"]["numerator_union"])
    m["SensCrownCalCOneBPct"] = macro_pct(cc["rungs"]["C1b"]["percent_union"])
    m["SensCrownCalNames"] = name_list(cc["partitions_removed"])
    m["SensCrownCalOwnNames"] = name_list(
        cc["partitions_removed_own_calibration"])
    m["SensCrownCalParentNames"] = name_list(
        cc["partitions_removed_parent_calibration"]) or "none"
    m["SensCrownUncalNames"] = name_list(cc["partitions_kept"])
    m["CrownCountedNoCalNames"] = name_list(bycal["none"])
    m["CrownCountedOneSidedNames"] = name_list(bycal["one_sided"])
    m["CrownCountedTwoSidedNames"] = name_list(bycal["two_sided"])
    # F3-12: the group that carries most of interval coverage, and its share
    c1b_by = sorted(((as_float(r["c1_units"]), r["partition"]) for r in basis
                     if is_euk(r) and r.get("interval_recoverable") == "yes"
                     and r.get("is_chronogram") == "yes"
                     and r.get("time_scale") == "absolute"), reverse=True)
    m["COneBLargestName"] = c1b_by[0][1]
    m["COneBLargestN"] = macro_int(c1b_by[0][0])
    m["COneBLargestSharePct"] = "%.0f" % (100.0 * c1b_by[0][0]
                                          / sum(x for x, _ in c1b_by))
    rp = s["c1b_rule_p_counted"]
    m["SensRulePCOneBN"] = macro_int(rp["rungs"]["C1b"]["numerator_union"])
    m["SensRulePCOneBPct"] = macro_pct(rp["rungs"]["C1b"]["percent_union"])
    m["SensRulePNames"] = name_list(rp["partitions_added"])
    m["RuleTipsOverNamedPct"] = "%.0f" % (100 * TABLE1_TIPS_OVER_NAMED)
    m["TabWhichTreeLargerNote"] = larger_tree_note()
    # F4-07: a larger tree built on a molecular source tree is described as
    # its source describes it (Data/larger_trees_not_counted.csv)
    for rs in larger_trees_not_counted().values():
        for r in rs:
            if r.get("n_molecular_source"):
                k = crown_macro_name(r["bibkey"])
                m["LargerN" + k] = macro_int(as_int(r["n"]))
                m["LargerNMol" + k] = macro_int(as_int(r["n_molecular_source"]))
    # every denominator lowered and raised by the elasticity fraction
    f_el = ELASTICITY_FRACTION
    m["RuleElasticityPct"] = "%.0f" % (100 * f_el)
    nob = s["leave_assembly_class_out"]["drop_automated_synthesis"]["rungs"]
    m["SensNoBirdsCZeroPct"] = macro_pct(nob["C0"]["percent_union"])
    m["SensNoBirdsCOnePct"] = macro_pct(nob["C1"]["percent_union"])

    ot = head["opentree_P01"]
    m["OTPhyloEukN"] = macro_int(ot["ot_phylo_extant_species_eukaryota"])
    m["OTPhyloEukPct"] = macro_pct(ot["ot_phylo_extant_species_eukaryota_pct"])
    m["OTCOneRatio"] = "%.2f" % ot["atlas_C1_over_ot_phylo"]
    m["OTKeysBoth"] = macro_int(ot["keys"]["keys_both"])
    me = ot["molecular_enumerated"]
    m["OTKeysMolecular"] = macro_int(me["keys_atlas_molecular"])
    m["OTKeysMolecularAndOT"] = macro_int(me["keys_molecular_and_ot"])
    m["OTAccMolecular"] = macro_int(me["accepted_molecular"])
    m["OTAccMolecularAndOT"] = macro_int(me["accepted_molecular_and_ot"])
    m["OTAccMolecularNotOT"] = macro_int(me["accepted_molecular_not_in_ot"])
    m["OTAccMolecularAndOTPct"] = macro_pct(pct(me["accepted_molecular_and_ot"],
                                                me["accepted_molecular"]))
    m["OTNPartitionsCompared"] = macro_int(len(me["partitions"]))
    # round F2 (F2-38): each overlap over its own base, in both directions
    m["OTKeysMolecularAndOTPct"] = macro_pct(pct(me["keys_molecular_and_ot"],
                                                 me["keys_atlas_molecular"]))
    m["OTKeysOTPhyloCompared"] = macro_int(me["keys_ot_phylo"])
    m["OTKeysMolecularAndOTOfOTPct"] = macro_pct(pct(me["keys_molecular_and_ot"],
                                                     me["keys_ot_phylo"]))
    # round F3 (F3-08): how the Open Tree species of the compared groups, of
    # the groups whose molecular tips are not listed, and of no group add up
    # to the eukaryotic count
    pp = ot["per_partition"]
    cmp_ = sum(pp[x]["ot_phylo_extant_species"] for x in me["partitions"])
    notc = sorted(x for x in pp if x not in me["partitions"])
    m["OTPhyloCompared"] = macro_int(cmp_)
    m["OTPhyloNotComparedN"] = macro_int(sum(pp[x]["ot_phylo_extant_species"]
                                             for x in notc))
    m["OTPhyloNotComparedNames"] = name_list(notc)
    m["OTPhyloOutsideGroups"] = macro_int(
        ot["ot_phylo_extant_species_eukaryota"]
        - ot["ot_phylo_extant_species_in_atlas_clades"])

    ca = head["crown_ages"]
    m["CrownNChron"] = macro_int(ca["n_absolute_chronograms"])
    m["CrownNCounted"] = macro_int(ca["n_with_counted_interval"])
    # the sets the text counts are also named, so that two sets of equal size
    # are never confused (round-3 finding)
    m["CrownNamesCounted"] = name_list(sorted(
        p for p, v in ca["per_partition"].items()
        if v["counts_for_c1b"] == "yes"))
    for k, v in ca["n_excluded_from_c1b"].items():
        m["CrownExcl" + crown_macro_name(k)] = macro_int(v)
    m["CrownNWithInterval"] = macro_int(sum(
        1 for v in ca["per_partition"].values() if v["interval_kind"]))
    m["CrownNNoInterval"] = macro_int(sum(
        1 for v in ca["per_partition"].values() if not v["interval_kind"]))
    # round F1 (F1-01, F1-07): which archives hold a tree set and which one
    # summary tree with an interval
    cset = sorted(p for p, v in ca["per_partition"].items()
                  if v["counts_for_c1b"] == "yes" and v["interval_kind"]
                  in CROWN_SET_KINDS)
    csum = sorted(p for p, v in ca["per_partition"].items()
                  if v["counts_for_c1b"] == "yes" and v["interval_kind"]
                  not in CROWN_SET_KINDS)
    oset = sorted(p for p, v in ca["per_partition"].items()
                  if v["counts_for_c1b"] != "yes" and v["interval_kind"]
                  in CROWN_SET_KINDS)
    m["CrownNCountedSet"] = macro_int(len(cset))
    m["CrownNamesCountedSet"] = name_list(cset)
    m["CrownNCountedSummary"] = macro_int(len(csum))
    m["CrownNamesCountedSummary"] = name_list(csum)
    m["CrownNSetNotCounted"] = macro_int(len(oset))
    m["CrownNamesSetNotCounted"] = name_list(oset)
    # F4-04: tree sets whose crown interval a calibration rule excludes; the
    # classification is of the crown only, and the set can still be
    # integrated over for the ages of other nodes
    m["CrownNamesSetCalib"] = name_list(sorted(
        p for p in oset if ca["per_partition"][p]["c1b_exclusion"]
        not in ("", None, "kind")))
    # F4-13: of those, the posterior sets (their node-age variation at other
    # nodes was not assessed)
    m["CrownNamesPostCalib"] = name_list(sorted(
        p for p in oset if ca["per_partition"][p]["c1b_exclusion"]
        not in ("", None, "kind")
        and ca["per_partition"][p]["interval_kind"] == "posterior_sample"))
    # round F2: the groups behind each loss of dated coverage, named
    pp = ca["per_partition"]
    m["CrownNamesNoInterval"] = name_list(sorted(
        p for p, v in pp.items() if not v["interval_kind"]))
    m["CrownNamesKind"] = name_list(sorted(
        p for p, v in pp.items() if v["interval_kind"]
        and v["c1b_exclusion"] == "kind"))
    m["CrownNamesCalib"] = name_list(sorted(
        p for p, v in pp.items() if v["interval_kind"]
        and v["c1b_exclusion"] not in ("", None, "kind")))

    # declared rule parameters the Methods state (rules, not results)
    import crown_intervals as _ci
    m["RuleChronCV"] = "%g" % CHRONOGRAM_TOLERANCE
    m["RuleDegenerateMyr"] = "%g" % _ci.DEGENERATE_WIDTH_MYR
    m["RuleBoundMyr"] = "%g" % _ci.BOUND_TOL_MYR
    m["RulePriorContainPct"] = "%g" % (100 * _ci.PRIOR_CONTAINMENT_MIN)
    m["RulePriorWidthRatio"] = "%g" % _ci.PRIOR_WIDTH_RATIO_MIN
    m["RuleIntervalPct"] = "%g" % round(100 * (_ci.QUANTILES[2]
                                                - _ci.QUANTILES[0]))

    # where C1 tips are lost on the way to C1b, by reason
    lost = coverage_losses(head, basis)
    for k, v in lost.items():
        m["Lost%sN" % k] = macro_int(v)
    m["LostUndatedOfCOnePct"] = macro_pct(pct(lost["NotChron"] + lost["RelTime"],
                                              rg["C1"]["numerator_union"]))
    m["LostRelTimeOfUndatedPct"] = macro_pct(pct(lost["RelTime"],
                                                 lost["NotChron"] + lost["RelTime"]))
    for k in ("NoInterval", "Kind", "Calib"):
        m["Lost%sOfCOneDPct" % k] = macro_pct(pct(lost[k],
                                                  rg["C1d"]["numerator_union"]))

    # the largest groups without a reference tree, named (R2-24)
    unrep = sorted((r for r in basis if is_euk(r)
                    and not as_int(r.get("c0_units"), 0)),
                   key=lambda r: (-as_int(r["described"], 0), r["partition"]))
    big = [r["partition"] for r in unrep[:UNREP_NAMED]]
    m["UnrepLargestNames"] = ", ".join(big[:-1]) + " and " + big[-1]

    # where the molecular species absent from Open Tree sit, by group (R2-10)
    ot = head["opentree_P01"]
    notot = ot["molecular_enumerated"]["accepted_molecular_not_in_ot"]
    for p in ot["molecular_enumerated"]["partitions"]:
        v = ot["per_partition"][p]
        m["OTNotOTSharePct" + crown_macro_name(p)] = macro_pct(
            pct(v["accepted_molecular_not_in_ot"], notot))

    # crown age over the root height of each fitted tree (R2-30, F3-04):
    # reported, not applied (D42)
    with open(WORKED_SEL, newline="", encoding="utf-8") as fh:
        sf = [float(r["crown_age_ma"]) / float(r["geometry_root_ma"])
              for r in csv.DictReader(fh)]
    m["WEScaleMin"] = "%.2f" % min(sf)
    m["WEScaleMax"] = "%.2f" % max(sf)

    w = head["worked_example_S13"]
    eff = w["rho_effect_same_tree"]
    fvm = w["full_vs_molecular"]
    fam = w["family_sampling_fit"]["per_partition"]
    m["WENEligible"] = macro_int(len(w["eligible"]))
    m["WENamesEligible"] = name_list(sorted(w["eligible"]))
    m["WENIdentity"] = macro_int(sum(v["classification"] == "identity"
                                     for v in eff.values()))
    m["WENBoundary"] = macro_int(sum(v["classification"] == "boundary"
                                     for v in eff.values()))
    m["WENRejected"] = macro_int(sum(bool(v["rho_full_rejected"])
                                     for v in eff.values()))
    m["WENFullGreater"] = macro_int(fvm["n_r_full_greater"])
    m["WENFullSmaller"] = macro_int(fvm["n_r_full_smaller"])
    m["WENFullOverlap"] = macro_int(fvm["n_r_intervals_overlap"])
    m["WENFamilyOverlap"] = macro_int(sum(bool(v["epsilon_intervals_overlap_uniform"])
                                          for v in fam.values()))
    # family against uniform sampling: the two likelihoods carry different
    # probabilities of sampling the same tips, so they are not compared; the
    # estimates are reported side by side (D30)
    m["WENFamilyEpsLower"] = macro_int(sum(
        v["epsilon"] < v["epsilon_uniform_tips"] for v in fam.values()))
    m["WENFamilyEpsHigher"] = macro_int(sum(
        v["epsilon"] > v["epsilon_uniform_tips"] for v in fam.values()))
    # F4-01 (D50): one test per eligible group, at the Bonferroni level over
    # the eligible groups, a family fixed before any fit is read; WENRejected
    # counts the summary trees rejected at that level
    tf = w["test_family"]
    m["RuleBonferroniLevel"] = "%.4f" % tf["level"]
    m["RuleBonferroniCritical"] = "%.2f" % tf["critical_value"]
    m["WENamesRejected"] = name_list(sorted(
        p for p, v in eff.items() if v["rho_full_rejected"])) or "none"
    m["RuleLRCritical"] = "%.2f" % WORKED_LR_CUTOFF
    m["RuleTestLevelPct"] = "%g" % (100 * WORKED_TEST_LEVEL)
    for p in w["eligible"]:
        n = crown_macro_name(p)
        e, f, fa = eff[p], fvm["per_partition"][p], fam[p]
        m["WERhoTips" + n] = "%.2f" % e["rho_tips"]
        m["WERhoFull" + n] = "%.2f" % e["rho_full"]
        m["WEEpsTips" + n] = "%.2f" % e["epsilon_rho_tips"]
        m["WEEpsFull" + n] = "%.2f" % e["epsilon_rho_full"]
        m["WELambdaRatio" + n] = "%.2f" % e["lambda_ratio_full_over_tips"]
        m["WELR" + n] = "%.2f" % e["lr_statistic"]
        # F3-15: the p-value of the statistic under the 50:50 mixture of a
        # point mass at zero and chi2(1) (Self & Liang 1987)
        if e["classification"] == "boundary":
            m["WEP" + n] = p_fmt(lr_mixture_p(e["lr_statistic"]))
        m["WERFull" + n] = "%.3f" % f["r_full"]
        m["WERMol" + n] = "%.3f" % f["r_molecular"]
        m["WERRatio" + n] = "%.2f" % f["r_ratio_full_over_molecular"]
        m["WEEpsFamily" + n] = "%.2f" % fa["epsilon"]
        # F4-12: the edge of the restriction, 1 - rho_m/rho_f, and the
        # profile interval of the matching extinction fraction on the
        # summary tree, which can include it
        m["WEEdge" + n] = "%.2f" % (1.0 - e["rho_tips"] / e["rho_full"])
        ftm = w["fits"][p + "|molecular|tips"]
        m["WEEpsTipsLo" + n] = "%.2f" % as_float(ftm["epsilon_lo"])
        m["WEEpsTipsHi" + n] = "%.2f" % as_float(ftm["epsilon_hi"])
        m["WERFamily" + n] = "%.3f" % fa["r"]
    # round F1 (F1-14): polytomies of the fitted trees and the fits with them
    # resolved at random
    po = w["polytomies"]
    with_poly = [p for p in w["eligible"]
                 if po[p]["events"]["full"]["n_polytomy_events"]
                 or po[p]["events"]["molecular"]["n_polytomy_events"]]
    m["WEPolyNGroups"] = macro_int(len(with_poly))
    m["WEPolyNames"] = name_list(sorted(with_poly))
    for p in with_poly:
        n = crown_macro_name(p)
        m["WEPolyPctFull" + n] = macro_pct(
            po[p]["events"]["full"]["polytomy_events_pct"])
        m["WEPolyPctPruned" + n] = macro_pct(
            po[p]["events"]["molecular"]["polytomy_events_pct"])
        lrs = [v["lr_statistic"] for v in po[p]["resolved"].values()]
        if eff[p]["classification"] == "boundary":
            m["WEPolyLRMin" + n] = "%.1f" % min(lrs)
            m["WEPolyLRMax" + n] = "%.1f" % max(lrs)
        rr = [v["r_ratio_full_over_molecular"] for v in po[p]["resolved"].values()]
        m["WEPolyRRatioMin" + n] = "%.2f" % min(rr)
        m["WEPolyRRatioMax" + n] = "%.2f" % max(rr)
    allr = [v["r_ratio_full_over_molecular"] for p in w["eligible"]
            for v in po[p]["resolved"].values()]
    m["WEPolyRRatioMin"] = "%.2f" % min(allr)
    m["WEPolyRRatioMax"] = "%.2f" % max(allr)
    m["WEPolyNBoundary"] = macro_int(sum(
        all(v["boundary"] for v in po[p]["resolved"].values())
        for p in w["eligible"]))
    m["WEPolyNRejected"] = macro_int(sum(
        all(v["rejected"] for v in po[p]["resolved"].values())
        for p in w["eligible"]))
    # round F2: the fits repeated over the archived tree sets
    ts = w["tree_sets"]
    for v in ts.values():
        n = crown_macro_name(v["partition"]) + WORKED_SET_TAG[v["set"]][0]
        m["WESetN" + n] = macro_int(v["n_trees"])
        m["WESetEpsMin" + n] = "%.2f" % v["epsilon_matching_min"]
        m["WESetEpsMax" + n] = "%.2f" % v["epsilon_matching_max"]
        m["WESetLRMin" + n] = "%.1f" % v["lr_min"]
        m["WESetLRMax" + n] = "%.1f" % v["lr_max"]
        m["WESetNBoundary" + n] = macro_int(v["n_boundary"])
        m["WESetLRMedian" + n] = "%.2f" % v["lr_median"]
        if v["r_ratio_min"] is not None:
            m["WESetRRatioMin" + n] = "%.2f" % v["r_ratio_min"]
            m["WESetRRatioMax" + n] = "%.2f" % v["r_ratio_max"]
    # F4-02, F4-03: for a group whose source archives a posterior set, the
    # distribution of the statistic over that set is the primary result; it
    # is described against the corrected critical value (where the statistic
    # lies, tree by tree), not summarised as a count of rejections
    crit = tf["critical_value"]
    post = {v["partition"]: v for v in ts.values()
            if v["set"] == WORKED_POSTERIOR_SET}
    m["WEPostNames"] = name_list(sorted(post))
    m["WESetNamesAboveAll"] = name_list(sorted(
        p for p, v in post.items() if v["lr_min"] > crit)) or "none"
    m["WESetNamesStraddle"] = name_list(sorted(
        p for p, v in post.items() if v["lr_min"] <= crit < v["lr_max"])) or "none"
    m["WESetNamesBelowAll"] = name_list(sorted(
        p for p, v in post.items() if v["lr_max"] <= crit)) or "none"
    nopost = [p for p in w["eligible"] if p not in post]
    m["WENoPostN"] = macro_int(len(nopost))
    m["WENoPostNames"] = name_list(sorted(nopost))
    m["WENoPostNRejected"] = macro_int(sum(
        bool(eff[p]["rho_full_rejected"]) for p in nopost))
    m["WESetNGroups"] = macro_int(len({v["partition"] for v in ts.values()}))
    m["WESetNamesGroups"] = name_list(sorted({v["partition"] for v in ts.values()}))
    allrr = [x for v in ts.values() if v["r_ratio_min"] is not None
             for x in (v["r_ratio_min"], v["r_ratio_max"])]
    m["WESetRRatioMin"] = "%.2f" % min(allrr)
    m["WESetRRatioMax"] = "%.2f" % max(allrr)

    for r in basis:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        n = crown_macro_name(r["partition"])
        d = as_int(r["described"], 0)
        c0 = as_int(r["c0_units"], 0)
        c1 = c0 * rung_weight(r, "C1")
        m["Described" + n] = macro_int(d)
        m["CZeroN" + n] = macro_int(c0)
        m["CZeroPct" + n] = macro_pct(pct(c0, d))
        m["COneN" + n] = macro_int(c1)
        m["COnePct" + n] = macro_pct(pct(c1, d))
        m["FMol" + n] = "%.2f" % as_float(r["f_mol"])

    out = ["% generated by scripts/rebuild_atlas.py from "
           "Paper/inputs/headline_numbers.json and Data/coverage_basis.csv.",
           "% Percentages: one decimal; below 1, two decimals.  Do not edit."]
    out += ["\\newcommand{\\%s}{%s}" % (k, v) for k, v in m.items()]
    return "\n".join(out) + "\n"


# ---------------------------------------------- revised manuscript inputs ---
#
# The revised MEE text (Paper/shared/manuscript_body_v2.tex) reads three
# generated tables and one figure projection.  Each builder reads only
# headline_numbers.json (as `head`), the basis table and the crown-age
# register, so the tables and the prose macros cannot disagree.

MOLECULAR_ONLY_TREES = os.path.join(ROOT, "Data", "molecular_only_trees.csv")
STUDY_DATING = OrderedDict([
    ("absolute_deposited", "absolute chronogram deposited"),
    ("absolute_deposited_fewer_tips", "dated; chronogram lacks tips"),
    ("absolute_not_deposited", "dated; chronogram not archived"),
    ("absolute_not_read", "dated; archive not readable"),
    ("relative_only", "relative time only"),
    ("not_dated", "not dated"),
])
STUDY_DATING_MACROS = OrderedDict([
    ("absolute_deposited_fewer_tips", "FewerTips"),
    ("absolute_not_deposited", "NotDeposited"),
    ("absolute_not_read", "NotRead"),
    ("relative_only", "Relative"),
    ("not_dated", "NotDated"),
])


def _molecular_only_trees():
    """R1-06: Data/molecular_only_trees.csv, the tree of molecular species each
    source distributes beside its taxonomically completed tree."""
    out = OrderedDict()
    for r in read_csv(MOLECULAR_ONLY_TREES):
        t = r["distributed_tree"].strip()
        if not t.startswith("https://") and not os.path.exists(
                os.path.join(ROOT, "Data", t)):
            die("%s: %s names %s, which is not under Data/"
                % (rel(MOLECULAR_ONLY_TREES), r["partition"], t))
        out[r["partition"]] = t
    return out


MOLONLY_BRANCH_CLASSES = ("time", "substitutions", "none")


def _molecular_only_dating():
    """Round F3 (F3-16).  What the tree of sequenced species each source
    distributes carries as branch lengths, measured on the file
    (Data/molecular_only_trees.csv columns branch_lengths, roottotip_cv,
    measured_file): 'none' when every branch length is zero or absent,
    'time' when the root-to-tip CV is within the chronogram rule
    (CHRONOGRAM_TOLERANCE), 'substitutions' otherwise.  A file held here is
    re-measured and must agree with the register; a file held elsewhere
    (the mammal DNA-only set) keeps its committed measurement."""
    out = OrderedDict()
    for r in read_csv(MOLECULAR_ONLY_TREES):
        p = r["partition"]
        klass = (r.get("branch_lengths") or "").strip()
        if klass not in MOLONLY_BRANCH_CLASSES:
            die("%s: %s has branch_lengths %r"
                % (rel(MOLECULAR_ONLY_TREES), p, klass))
        f = os.path.join(ROOT, "Data", (r.get("measured_file") or "").strip())
        if r.get("measured_file") and os.path.exists(f):
            text = read_first_tree(f)
            root = PARSE_NEWICK(normalize_newick(text))
            stack, nonzero = [(root, True)], 0
            while stack:
                (kids, _, bl), top = stack.pop()
                if not top and bl == bl and bl != 0:
                    nonzero += 1
                for k in kids:
                    stack.append((k, False))
            _, depths = tree_facts(text)
            cv = coefficient_of_variation(depths) if nonzero else None
            got = ("none" if not nonzero
                   else "time" if cv <= CHRONOGRAM_TOLERANCE
                   else "substitutions")
            if got != klass or (cv is not None and abs(
                    cv - as_float(r["roottotip_cv"], -1)) > CV_RECHECK_TOL):
                die("%s: %s measures %s (CV %s) on %s; the register says %s "
                    "(CV %s)" % (rel(MOLECULAR_ONLY_TREES), p, got, cv,
                                 r["measured_file"], klass, r["roottotip_cv"]))
        out[p] = klass
    return out


def _study_dating():
    """R1-09: what the study behind each tree did about dating
    (Data/dating_deposit_checks.csv, column study_dating)."""
    out = {}
    for r in read_csv(DATING_DEPOSIT_CHECKS):
        v = (r.get("study_dating") or "").strip()
        if v not in STUDY_DATING:
            die("%s: %s has study_dating %r" % (rel(DATING_DEPOSIT_CHECKS),
                                                 r["partition"], v))
        out[r["partition"]] = v
    return out


WHICH_TREE_HEADER = ["partition", "category", "bibkey", "flags", "tips",
                     "named_pct", "molecular_tips", "molecular_basis",
                     "molecular_pct", "time_scale", "crown_age_ma",
                     "interval_lo_ma", "interval_hi_ma", "interval_kind",
                     "interval_status", "rho_full", "rho_molecular",
                     "described", "study_dating", "molecular_only_tree",
                     "molecular_only_dated",
                     "molecular_share", "tips_exceed_named",
                     "larger_tree_not_counted"]
# Round F1 (F1-11, F1-17): Table 1 flags a group whose tips exceed its named
# species by more than this fraction, where tips over described species
# overstates the share of the group's species that the tree samples.
TABLE1_TIPS_OVER_NAMED = 0.10
INTERVAL_KIND_TEXT = OrderedDict([
    ("posterior_sample", "posterior tree set"),
    ("node_hpd_annotation", "HPD on summary tree"),
    ("node_label_interval", "intervals on tree"),
    ("bootstrap_dating_sample", "bootstrap dated set"),
    ("calibration_prior_sample", "set under calibration prior"),
    ("taxon_addition_sample", "taxon-addition set"),
    ("date_assignment_sample", "date-assignment set"),
])
INTERVAL_STATUS_TEXT = OrderedDict([
    ("counted", "counted"),
    ("kind", "other variation"),
    ("degenerate", "zero width"),
    ("at_bound", "at calibration bound"),
    ("prior_reproduced", "within prior, not narrowed"),
    ("prior_sampled", "root ages drawn from prior"),
    ("none", "none archived"),
])
TIME_SCALE_TEXT = {"absolute": "absolute", "relative": "relative",
                   "not_ultrametric": "not a chronogram"}
# Table 1 prints one- or two-letter codes so that no cell wraps; the caption
# (Paper/shared/metadata_v2.tex) decodes every code and verify_atlas.py checks
# that each code printed is decoded there (round-3 layout finding).
TABLE1_DATES_CODE = OrderedDict([
    ("absolute", "A"), ("relative", "R"),
    ("relative_only", "S"),
    ("absolute_deposited_fewer_tips", "Dt"),
    ("absolute_not_deposited", "Dn"),
    ("absolute_not_read", "Du"),
    ("not_dated", "N"),
])
TABLE1_INTERVAL_CODE = OrderedDict([
    ("counted", "u"), ("kind", "v"), ("degenerate", "z"), ("at_bound", "b"),
    ("prior_reproduced", "p"), ("prior_sampled", "r"), ("none", "n"),
])


def name_list(names):
    """'A, B and C': group names in running text."""
    names = list(names)
    if not names:
        return ""
    return names[0] if len(names) == 1 else \
        ", ".join(names[:-1]) + " and " + names[-1]


def rho_fmt(x):
    """A sampling fraction: two decimals, or two significant digits below
    0.1 so that a small group's fraction does not print as 0.00."""
    return "%.2f" % x if x >= 0.1 else "%#.2g" % x


def _molecular_tip_counts():
    n = Counter()
    for r in read_csv(os.path.join(ROOT, "Data", "molecular_tips.csv")):
        n[r["partition"]] += 1
    return n


def _which_tree_rows(head, basis):
    reg = {r["filename"]: r for r in read_csv(DATASET_REGISTRY)}
    mol = _molecular_tip_counts()
    crown = head["crown_ages"]["per_partition"]
    sd = _study_dating()
    molonly = _molecular_only_trees()
    modat = _molecular_only_dating()
    larger = larger_trees_not_counted()
    rows = []
    for r in basis:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        p = r["partition"]
        g = reg.get(r["canonical_tree"])
        if g is None:
            die("Table 1: %s is not in %s" % (r["canonical_tree"],
                                               rel(DATASET_REGISTRY)))
        d = as_int(r["described"], 0)
        tips = as_int(r["labels_unique"], 0)
        if r["c1_basis"] == "enumerated":
            mtips = mol[p]
        elif r["c1_basis"] == "all_tips":
            mtips = tips
        else:
            die("Table 1: %s has c1_basis %r" % (p, r["c1_basis"]))
        c = crown.get(p)
        if r["is_chronogram"] == "yes" and c is None:
            die("Table 1: chronogram %s has no crown-age register row" % p)
        kind = c["interval_kind"] if c else ""
        if not c or not kind:
            status = "none"
        elif c["counts_for_c1b"] == "yes":
            status = "counted"
        else:
            status = c["c1b_exclusion"]
        rows.append(OrderedDict([
            ("partition", p), ("category", r["category"]),
            ("bibkey", g["bibkey"]), ("flags", g.get("flags") or ""),
            ("tips", tips),
            ("named_pct", macro_pct(pct(as_int(r["c0_units"], 0), d))),
            ("molecular_tips", mtips), ("molecular_basis", r["c1_basis"]),
            ("molecular_pct", macro_pct(pct(as_float(r["c1_units"]), d))),
            ("time_scale", r["time_scale"]),
            ("crown_age_ma", "%.1f" % as_float(c["crown_age_ma"]) if c else ""),
            ("interval_lo_ma", "%.1f" % as_float(c["interval_lo_ma"])
             if kind else ""),
            ("interval_hi_ma", "%.1f" % as_float(c["interval_hi_ma"])
             if kind else ""),
            ("interval_kind", kind), ("interval_status", status
                                      if r["is_chronogram"] == "yes" else ""),
            ("rho_full", rho_fmt(tips / d)),
            ("rho_molecular", rho_fmt(mtips / d)),
            ("described", d),
            ("study_dating", sd.get(p) or (
                "absolute_deposited" if r["is_chronogram"] == "yes"
                and r["time_scale"] == "absolute" else "relative_only"
                if r["time_scale"] == "relative" else "")),
            ("molecular_only_tree", "yes" if p in molonly else "no"),
            # F3-16: whether that tree is dated (measured on the file)
            ("molecular_only_dated", ("yes" if modat.get(p) == "time"
                                      else "no") if p in molonly else ""),
            # molecular tips over tips: the molecular fraction of the tree,
            # which the text quotes (R2-17), beside the sampling fractions
            ("molecular_share", "%.2f" % (mtips / tips)),
            ("tips_exceed_named", "yes" if tips > (1 + TABLE1_TIPS_OVER_NAMED)
             * as_int(r["c0_units"], 0) else "no"),
            ("larger_tree_not_counted", ";".join(
                x["bibkey"] for x in larger.get(p, []))),
        ]))
        if not rows[-1]["study_dating"]:
            die("Table 1: %s is not a chronogram and has no study_dating "
                "row in %s" % (p, rel(DATING_DEPOSIT_CHECKS)))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             -float(x["named_pct"]), x["partition"]))
    return rows


def build_which_tree_csv(head, basis, warnings):
    return render_csv(WHICH_TREE_HEADER, _which_tree_rows(head, basis))


def build_which_tree_tex(head, basis, warnings):
    """Table 1 body rows (12 columns); the environment and caption are in the
    manuscript.  Superscripts: a molecular tips stated by the source, b an
    automated synthesis, c a preprint source, d the source also distributes a
    tree of its molecular species (R1-06).  The dates column says what the
    tree is and, for a tree that is not an absolute chronogram, what the
    study did about dating (R1-09); the interval column says whether the
    deposit's interval counts and, if not, why (crown ages are in Fig. 4)."""
    out, prev = [], None
    for r in _which_tree_rows(head, basis):
        if r["category"] != prev:
            if prev is not None:
                out.append(r"\midrule")
            out.append(r"\multicolumn{12}{@{}l}{\textit{%s}} \\"
                       % tex_escape(r["category"]))
            prev = r["category"]
        src = r"\citet{%s}" % r["bibkey"]
        if "preprint" in r["flags"].split(";"):
            src += r"$^{c}$"
        b = next(x for x in basis if x["partition"] == r["partition"])
        if b.get("topology_assembly") == ASSEMBLY_AUTOMATED:
            src += r"$^{b}$"
        if r["larger_tree_not_counted"]:
            src += r"$^{f}$"
        mt = thousands(r["molecular_tips"])
        if r["molecular_basis"] == "all_tips":
            mt += r"$^{a}$"
        if r["molecular_only_tree"] == "yes":
            mt += (r"$^{d}$" if r["molecular_only_dated"] == "yes"
                   else r"$^{g}$")
        if r["time_scale"] in ("absolute", "relative"):
            dates = TABLE1_DATES_CODE[r["time_scale"]]
        else:
            dates = TABLE1_DATES_CODE[r["study_dating"]]
        if r["interval_status"]:
            unc = TABLE1_INTERVAL_CODE[r["interval_status"]]
        else:
            unc = "---"
        rm = (r["rho_molecular"] if r["molecular_tips"] != r["tips"]
              else "=")
        rf = r["rho_full"] + (r"$^{e}$" if r["tips_exceed_named"] == "yes"
                              else "")
        out.append(" & ".join([
            tex_escape(r["partition"]), src, thousands(r["described"]),
            thousands(r["tips"]), r["named_pct"], mt, r["molecular_share"],
            r["molecular_pct"], dates, unc, rf, rm]) + r" \\")
    return "\n".join(out) + "\n"


SENS_COLUMNS = ["rep", "col", "union", "gbif", "product", "calib", "allint",
               "nobirds"]


def _rung_sensitivity_rows(head):
    names = (("C0", "named"), ("C1", "molecular"), ("C1d", "dated"),
             ("C1b", "interval"))
    rg, s = head["rungs"], head["sensitivities"]
    rep = s["denominator"]["represented_partitions_only"]["rungs"]
    col = s["denominator"]["col_single_authority"]["filled"]["rungs"]
    u = head["union_bound_S09"]
    prod = s["c1_definition"]["product"]
    cal = s["c1b_calibration_fixed"]["rungs"]
    alli = s["c1b_all_intervals"]["rungs"]
    nob = s["leave_assembly_class_out"]["drop_automated_synthesis"]["rungs"]

    def pu(e):
        return macro_pct(e.get("percent_union", e["percent"]))
    rows = []
    for key, label in names:
        rows.append(OrderedDict([
            ("level", label),
            ("species", macro_int(rg[key]["numerator_union"])),
            ("pct", macro_pct(rg[key]["percent_union"])),
            ("rep", pu(rep[key])),
            ("col", pu(col[key])),
            ("union", macro_pct(u["C0_union_pct"]) if key == "C0" else
             ("<= " + macro_pct(u["C1_union_high_pct"])) if key == "C1"
             else ""),
            ("gbif", macro_pct(s["c0_gbif_unmatched"]["percent"])
             if key == "C0" else ""),
            ("product", macro_pct(prod[key]["percent"])
             if key in ("C1", "C1d") else ""),
            ("calib", pu(cal["C1b"]) if key == "C1b" else ""),
            ("allint", pu(alli["C1b"]) if key == "C1b" else ""),
            ("nobirds", pu(nob[key])),
        ]))
    return rows


def build_rungs_sensitivity_csv(head, warnings):
    return render_csv(["level", "species", "pct"] + SENS_COLUMNS,
                      _rung_sensitivity_rows(head))


def build_rungs_sensitivity_tex(head, warnings):
    out = []
    for r in _rung_sensitivity_rows(head):
        cells = [r["level"].capitalize(), r["species"], r["pct"]]
        for c in SENS_COLUMNS:
            v = r[c]
            cells.append("---" if not v else
                         r"$\leq$\,%s" % v[3:] if v.startswith("<= ") else v)
        out.append(" & ".join(cells) + r" \\")
    return "\n".join(out) + "\n"


WORKED_HEADER = ["partition", "tips_full", "tips_molecular", "rho_molecular",
                 "rho_fulltree", "lambda_ratio", "eps_molecular",
                 "eps_molecular_lo", "eps_molecular_hi", "eps_fullrho",
                 "eps_fullrho_lo", "eps_fullrho_hi", "lr", "classification",
                 "r_full", "r_full_lo", "r_full_hi", "r_molecular",
                 "r_molecular_lo", "r_molecular_hi", "eps_family",
                 "eps_family_lo", "eps_family_hi", "r_family", "r_family_lo",
                 "r_family_hi", "basis", "n_trees", "n_boundary", "lr_lo",
                 "lr_hi"]
# Round F4 (F4-10, D54).  Where the source archives the dated trees its
# summary tree summarises, Table 3 and Figure 5 report the median and range
# over that set: the posterior set, else the set with dates drawn at random
# (Birds).  Consensus summaries carry polytomies whose ages come from the
# summary, and their profile intervals describe that one tree.  A set dated
# under root ages drawn from the prior (Sharks) varies the prior, not the
# data, and is not used.  The summary-tree fits stay in the Supplement.
WORKED_TABLE_SETS = ("posterior", "date_random")


def worked_table_set(p, rows=None):
    """(set name, rows of worked_example_bd_sets.csv) Table 3 uses for p,
    or (None, [])."""
    rows = read_csv(WORKED_SETS) if rows is None else rows
    for st in WORKED_TABLE_SETS:
        rs = [r for r in rows if r["partition"] == p and r["set"] == st]
        if rs:
            return st, rs
    return None, []


def _set_lambda_ratio(r):
    lm = as_float(r["r_matching"]) / (1.0 - as_float(r["epsilon_matching"]))
    lf = as_float(r["r_full_rho"]) / (1.0 - as_float(r["epsilon_full_rho"]))
    return lf / lm


def _worked_rows(head):
    w = head["worked_example_S13"]
    fits, eff = w["fits"], w["rho_effect_same_tree"]
    srows = read_csv(WORKED_SETS)
    rows = []
    for p in w["eligible"]:
        ft, fc = fits[p + "|molecular|tips"], fits[p + "|molecular|full_misapplied"]
        ff, fm = fits[p + "|full|tips"], fits[p + "|molecular|family"]

        def f2(x):
            return "%.2f" % as_float(x)

        def f3(x):
            return "%.3f" % as_float(x)
        rows.append(OrderedDict([
            ("partition", p), ("tips_full", as_int(ff["n_tips"])),
            ("tips_molecular", as_int(ft["n_tips"])),
            ("rho_molecular", f2(ft["rho"])), ("rho_fulltree", f2(fc["rho"])),
            ("lambda_ratio", "%.2f" % eff[p]["lambda_ratio_full_over_tips"]),
            ("eps_molecular", f2(ft["epsilon"])),
            ("eps_molecular_lo", f2(ft["epsilon_lo"])),
            ("eps_molecular_hi", f2(ft["epsilon_hi"])),
            ("eps_fullrho", f2(fc["epsilon"])),
            ("eps_fullrho_lo", f2(fc["epsilon_lo"])),
            ("eps_fullrho_hi", f2(fc["epsilon_hi"])),
            ("lr", "%.2f" % eff[p]["lr_statistic"]),
            ("classification", eff[p]["classification"]
             + (" (rejected)" if eff[p]["rho_full_rejected"] else "")),
            ("r_full", f3(ff["r"])), ("r_full_lo", f3(ff["r_lo"])),
            ("r_full_hi", f3(ff["r_hi"])),
            ("r_molecular", f3(ft["r"])), ("r_molecular_lo", f3(ft["r_lo"])),
            ("r_molecular_hi", f3(ft["r_hi"])),
            ("eps_family", f2(fm["epsilon"])),
            ("eps_family_lo", f2(fm["epsilon_lo"])),
            ("eps_family_hi", f2(fm["epsilon_hi"])),
            ("r_family", f3(fm["r"])), ("r_family_lo", f3(fm["r_lo"])),
            ("r_family_hi", f3(fm["r_hi"])),
            ("basis", "summary"), ("n_trees", ""), ("n_boundary", ""),
            ("lr_lo", ""), ("lr_hi", ""),
        ]))
        st, rs = worked_table_set(p, srows)
        if st:
            row = rows[-1]

            def put(key, col, fmt):
                xs = [as_float(r[col]) for r in rs]
                row[key] = fmt(median(xs))
                row[key + "_lo"], row[key + "_hi"] = fmt(min(xs)), fmt(max(xs))
            put("eps_molecular", "epsilon_matching", f2)
            put("eps_fullrho", "epsilon_full_rho", f2)
            put("r_full", "r_full_tree", f3)
            put("r_molecular", "r_matching", f3)
            put("lr", "lr_statistic", lambda x: "%.2f" % x)
            nb = sum(1 for r in rs
                     if as_float(r["epsilon_full_rho"]) < WORKED_EPS_BOUNDARY
                     and as_float(r["lr_statistic"]) > 2 * WORKED_IDENTITY_LL_TOL)
            row["lambda_ratio"] = "%.2f" % median([_set_lambda_ratio(r) for r in rs])
            row["classification"] = ("boundary" if nb == len(rs) else
                                     "identity" if nb == 0 else "mixed")
            row["basis"], row["n_trees"], row["n_boundary"] = st, str(len(rs)), str(nb)
    return rows


def build_worked_table_csv(head, warnings):
    return render_csv(WORKED_HEADER, _worked_rows(head))


def build_worked_table_tex(head, warnings):
    """Table 3 body: one row per eligible group, 13 columns.  The pair of
    fits is printed as 'rescaled' (the identity holds) or 'mu = 0' (the
    misapplied fraction forces extinction to zero)."""
    out = []
    mark = {"posterior": r"$^{\dagger}$", "date_random": r"$^{\ddagger}$"}
    for r in _worked_rows(head):
        def ci(k):
            return r"%s (%s--%s)" % (r[k], r[k + "_lo"], r[k + "_hi"])
        if r["basis"] == "summary":
            cls = r["classification"].replace(" (rejected)", r"$^{*}$")
            cls = cls.replace("identity", "rescaled").replace(
                "boundary", r"$\mu = 0$")
            lr = r["lr"] + " " + cls
            fam = (ci("eps_family"), ci("r_family"))
        else:
            # F4-10, F4-11: over the set, the statistic's median (range) and
            # the trees in which the restriction binds; the family fit is of
            # the summary tree only ($^{s}$)
            lr = r"%s (%s--%s), $\mu = 0$ in %s/%s" % (
                r["lr"], r["lr_lo"], r["lr_hi"], r["n_boundary"], r["n_trees"])
            fam = (ci("eps_family") + r"$^{s}$", ci("r_family") + r"$^{s}$")
        out.append(" & ".join([
            tex_escape(r["partition"]) + mark.get(r["basis"], ""),
            thousands(r["tips_full"]),
            thousands(r["tips_molecular"]), r["rho_molecular"],
            r["rho_fulltree"], r["lambda_ratio"], ci("eps_molecular"),
            ci("eps_fullrho"), lr, ci("r_full"),
            ci("r_molecular"), fam[0], fam[1]]) + r" \\")
    return "\n".join(out) + "\n"


CASCADE_HEADER = ["kind", "level", "group", "species"]


LEVEL_ORDER = ("C0", "C1", "C1d", "C1b")
UNREP_NAMED = 4     # how many of the largest groups without a tree are named
LOSS_REASONS = ("NotChron", "RelTime", "NoInterval", "Kind", "Calib")


def union_level_assignment(head, basis, shipped):
    """The four levels and the losses between them on the union convention
    (round-2 reader R2-15): every accepted species is counted once.  A
    species named in the trees of several groups is assigned to the group in
    which it reaches the furthest level (ties: the group first in name
    order), so the category stacks of Figure 1a sum to the union totals and
    the losses of Figure 1b are set differences of consecutive unions.  The
    reason a species is lost is that of its assigned group: NotChron or
    RelTime (the tree is not a chronogram in absolute time), NoInterval,
    Kind (the archived set describes other variation) or Calib (the interval
    is set by a calibration)."""
    tr, ca = head["trees"], head["crown_ages"]["per_partition"]
    reltime = set(tr["relative_time_partitions"])
    groups = {}
    for r in basis:
        p = r["partition"]
        if not is_euk(r) or not as_int(r.get("c0_units"), 0):
            continue
        c1 = shipped.c1_names.get(p)
        if c1 is None or float(len(c1)) != as_float(r.get("c1_units")):
            die("%s: the union losses need an enumerated C1 set" % p)
        chron = r.get("is_chronogram") == "yes"
        v = ca.get(p) or {}
        b = chron and v.get("counts_for_c1b") == "yes"
        if not chron:
            reason = "RelTime" if p in reltime else "NotChron"
        elif b:
            reason = ""
        elif not v.get("interval_kind"):
            reason = "NoInterval"
        elif v.get("c1b_exclusion") == "kind":
            reason = "Kind"
        else:
            reason = "Calib"
        groups[p] = (r["category"], set(shipped.c0_names[p]), set(c1),
                     chron, b, reason)
    best = {}
    for p in sorted(groups):
        cat, c0, c1, chron, b, reason = groups[p]
        for sp in c0:
            lvl = (3 if (sp in c1 and b) else 2 if (sp in c1 and chron)
                   else 1 if sp in c1 else 0)
            if sp not in best or lvl > best[sp][0]:
                best[sp] = (lvl, p)
    levels = OrderedDict((k, OrderedDict((c, 0) for c in CATEGORY_ORDER))
                         for k in LEVEL_ORDER)
    lost = OrderedDict([("TaxPlaced", 0)] + [(k, 0) for k in LOSS_REASONS])
    shared = 0
    for sp, (lvl, p) in best.items():
        cat, reason = groups[p][0], groups[p][5]
        for i in range(lvl + 1):
            levels[LEVEL_ORDER[i]][cat] += 1
        if lvl == 0:
            lost["TaxPlaced"] += 1
        elif lvl < 3:
            lost[reason] += 1
    counts = Counter(sp for g in groups.values() for sp in g[1])
    shared = sum(1 for c in counts.values() if c > 1)
    return OrderedDict([
        ("rule", "union convention: each accepted species once, assigned to "
                 "the group in which it reaches the furthest level (ties: "
                 "first group in name order); losses are set differences of "
                 "consecutive levels, by the reason of the assigned group"),
        ("species_in_several_groups", shared),
        ("levels_by_category", levels),
        ("lost", lost),
    ])


def coverage_losses(head, basis):
    """Molecular coverage lost on the way to interval coverage, by reason, on
    the union convention (union_level_assignment)."""
    lost = head["coverage_union"]["lost"]
    return OrderedDict((k, lost[k]) for k in LOSS_REASONS)


def build_fig_cascade_csv(head, basis, warnings):
    """Figure 1 input: the four levels by category, and the losses between
    them by reason, on the union convention (eukaryotic only)."""
    cu = head["coverage_union"]
    cats = {x["category"] for x in basis if is_euk(x)}
    rows = []
    for lvl in LEVEL_ORDER:
        for c, v in cu["levels_by_category"][lvl].items():
            if c in cats:
                rows.append(OrderedDict([("kind", "level"), ("level", lvl),
                                         ("group", c), ("species", "%d" % v)]))
    for k, v in cu["lost"].items():
        step = ("C0>C1" if k == "TaxPlaced" else
                "C1>C1d" if k in ("NotChron", "RelTime") else "C1d>C1b")
        rows.append(OrderedDict([("kind", "loss"), ("level", step),
                                 ("group", k), ("species", "%d" % v)]))
    return render_csv(CASCADE_HEADER, rows)


def union_bound(truth, basis, shipped, warnings):
    """S09.  C0 and C1 over the UNION of all distributed trees, a second bound
    on the canonical-tree figures.

    A tip of any shipped eukaryotic tree enters C0 when it resolves to a
    Catalogue of Life accepted name (the same rule as the canonical tips,
    Data/standardized/col_resolution_noncanonical.csv).  It is assigned to a
    partition AS A TIP: to the partition whose canonical tips carry its CoL
    family most often (the family map, measured from col_resolution.csv);
    when its family occurs in no canonical tree, to its tree's partition,
    which is truth_dated.csv's partition_name or, for a tree tagged by its
    compiler, the majority partition of its mapped tips.  An accepted name
    already counted by any canonical tree is not counted again.

    C1 over the union is bounded, not estimated: the molecular status of a
    tip outside the canonical trees is not measured, so the added names enter
    with weight 0 (low) and 1 (high)."""
    canon_parts = {r["partition"] for r in basis
                   if is_euk(r) and as_int(r["c0_units"], 0) > 0}
    fam = defaultdict(Counter)
    for part in canon_parts:
        for r in shipped.col_rows.get(part, {}).values():
            if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_family"]:
                fam[r["col_family"]][part] += 1
    fam_part = {f: sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                for f, c in fam.items()}
    counted = set()
    for part in canon_parts:
        counted |= set(shipped.c0_names.get(part, {}))
    rows = defaultdict(list)
    if os.path.exists(COL_RESOLUTION_NONCANON):
        with open(COL_RESOLUTION_NONCANON, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                rows[r["filename"]].append(r)
    else:
        die("%s is missing; run scripts/match_col.py"
            % rel(COL_RESOLUTION_NONCANON))
    noncanon = [fn for fn, d in truth.dated.items()
                if not truthy(d.get("is_canonical"))
                and not ((d.get("partition_name") or "").strip()
                         in truth.species
                         and truth.is_prokaryotic(d["partition_name"].strip()))]
    missing = sorted(set(noncanon) - set(rows))
    if missing:
        die("non-canonical tree(s) with no row in %s: %s"
            % (rel(COL_RESOLUTION_NONCANON), ", ".join(missing[:5])))
    added = defaultdict(dict)       # partition -> {col_id: first file}
    per_file = OrderedDict()
    for fn in sorted(noncanon):
        res = [r for r in rows[fn] if r["match_type"] in COL_RESOLUTION_RESOLVED
               and r["col_id"]]
        tag = (truth.dated[fn].get("partition_name") or "").strip()
        if tag in canon_parts:
            tree_part = tag
        else:
            c = Counter(fam_part[r["col_family"]] for r in res
                        if r["col_family"] in fam_part)
            tree_part = (sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                         if c else "")
        n_new = 0
        for r in res:
            part = fam_part.get(r["col_family"]) or tree_part
            if not part or r["col_id"] in counted:
                continue
            if r["col_id"] not in added[part]:
                added[part][r["col_id"]] = fn
                n_new += 1
        per_file[fn] = OrderedDict([("tree_partition", tree_part),
                                    ("resolved_tips", len(res)),
                                    ("names_added", n_new)])
    union_ids = set().union(*added.values()) if added else set()
    den = sum(as_int(r["described"], 0) for r in basis if is_euk(r))
    c0 = sum(as_int(r["c0_units"], 0) for r in basis if is_euk(r))
    c1 = round(sum(as_int(r["c0_units"], 0) * rung_weight(r, "C1")
                   for r in basis if is_euk(r)), NUMERATOR_DP)
    per_part = OrderedDict()
    for r in basis:
        part = r["partition"]
        if not is_euk(r):
            continue
        n0 = as_int(r["c0_units"], 0)
        add = len(added.get(part, {}))
        if not n0 and not add:
            continue
        den_p = as_int(r["described"], 0)
        per_part[part] = OrderedDict([
            ("C0_canonical", n0), ("names_added_by_other_trees", add),
            ("C0_union", n0 + add),
            ("C0_canonical_pct", pct(n0, den_p)),
            ("C0_union_pct", pct(n0 + add, den_p)),
        ])
    n_add = len(union_ids)
    return OrderedDict([
        ("definition",
         "C0 over the union of all distributed eukaryotic trees: accepted "
         "names of the canonical tips plus those of every other shipped "
         "tree's tips not already counted; each added tip assigned to a "
         "partition by its CoL family (map measured on the canonical tips), "
         "else by its tree's partition. C1 over the union is a band: added "
         "names at weight 0 and at weight 1."),
        ("n_noncanonical_trees", len(noncanon)),
        ("names_added", n_add),
        ("C0_canonical", c0),
        ("C0_union", c0 + n_add),
        ("C0_canonical_pct", pct(c0, den)),
        ("C0_union_pct", pct(c0 + n_add, den)),
        ("C1_canonical", c1),
        ("C1_union_low", c1),
        ("C1_union_high", round(c1 + n_add, NUMERATOR_DP)),
        ("C1_union_low_pct", pct(c1, den)),
        ("C1_union_high_pct", pct(c1 + n_add, den)),
        ("denominator", den),
        ("per_partition", per_part),
        ("per_file", per_file),
    ])


OPENTREE_COUNTS = os.path.join(ROOT, "Data", "external",
                               "opentree16.1_counts.csv")
OPENTREE_SPECIES = os.path.join(ROOT, "Data", "external",
                                "opentree16.1_phylo_species.tsv.gz")
# R1-14: the same Open Tree species resolved to COL26.9 accepted names by
# scripts/match_col.py --opentree (the rule of the atlas tips).
OPENTREE_SPECIES_COL = os.path.join(ROOT, "Data", "external",
                                    "opentree16.1_phylo_species_col.csv")
OT_HEADER = ["partition", "described", "atlas_C0", "atlas_C1",
             "atlas_C1_tips_enumerated", "ot_tips_extant_species",
             "ot_phylo_extant_species", "atlas_C1_over_ot_phylo",
             "keys_atlas", "keys_ot_phylo", "keys_both", "keys_atlas_only",
             "keys_ot_only", "keys_union", "keys_atlas_molecular",
             "keys_molecular_and_ot", "keys_molecular_union_ot",
             "accepted_molecular", "accepted_molecular_and_ot",
             "accepted_molecular_not_in_ot"]


def opentree_comparison(basis, shipped, mol):
    """P01.  The atlas against Open Tree of Life synthesis v16.1.

    Open Tree's counterpart of C1 is its 'tips from phylogeny': leaves of the
    grafted solution, i.e. taxa placed by at least one input phylogeny.  Only
    those that are extant species-rank taxa are compared (Data/external/,
    written by scripts/opentree_extract.py).  Open Tree counts a SET of
    leaves, including ones placed by morphological input trees; atlas C1 is an
    f_mol-weighted count.  Names are compared on the S05 query of each atlas
    tip (its species key, decision D4, or its documented override) against
    the Open Tree name's Genus_epithet, by exact spelling and without
    synonymy, so 'keys_both' is a lower bound."""
    for path in (OPENTREE_COUNTS, OPENTREE_SPECIES):
        if not os.path.exists(path):
            die("%s is missing; run scripts/opentree_extract.py" % rel(path))
    if not os.path.exists(OPENTREE_SPECIES_COL):
        die("%s is missing; run scripts/match_col.py --opentree"
            % rel(OPENTREE_SPECIES_COL))
    mol_sids = defaultdict(set)
    for r in read_csv(os.path.join(ROOT, "Data", "molecular_tips.csv")):
        mol_sids[r["partition"]].add(r["standardized_id"])
    ot_acc = defaultdict(set)
    for r in read_csv(OPENTREE_SPECIES_COL):
        if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
            ot_acc[r["clade"]].add(r["col_id"])
    with open(OPENTREE_COUNTS, newline="", encoding="utf-8") as fh:
        cnt = {r["clade"]: r for r in csv.DictReader(fh)}
    ot_keys = defaultdict(set)
    with gzip.open(OPENTREE_SPECIES, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            _u, nm, clade = line.rstrip("\n").split("\t")
            k, ok = SPECIES_KEY(nm.replace(" ", "_"), "ot")
            if ok:
                ot_keys[clade].add(k)
    rows = []
    for r in basis:
        part = r["partition"]
        if not is_euk(r) or part not in cnt or as_int(r["c0_units"], 0) == 0:
            continue
        idmap = shipped.ids.get(part, {})
        colq = shipped.col_rows.get(part, {})

        def key_of(nm):
            # the S05 query: the species key, or the documented override for
            # a label that carries no binomial (crustacean and crocodilian
            # codes); a non-species label (D13) has none
            sid = idmap.get(nm, nm)
            q = (colq.get(sid) or {}).get("query")
            if q is not None:
                return q or None
            k, ok = SPECIES_KEY(nm, sid)
            return k if ok else None

        akeys = {k for k in map(key_of, shipped.labels.get(part, ())) if k}
        ev = mol.get(part)
        mkeys = None
        macc = None
        if ev is not None and ev.tip_names is not None and ev.source_file:
            mkeys = {k for k in map(key_of, ev.tip_names) if k}
            # accepted names of the molecular tips of Data/molecular_tips.csv
            # (the C1 units of D19)
            macc = set()
            for sid in mol_sids.get(part, ()):
                cr = colq.get(sid) or {}
                if cr.get("match_type") in COL_RESOLUTION_RESOLVED \
                        and cr.get("col_id"):
                    macc.add(cr["col_id"])
        o = ot_keys.get(part, set())
        oa = ot_acc.get(part, set())
        c1 = round(as_int(r["c0_units"], 0) * rung_weight(r, "C1"),
                   NUMERATOR_DP)
        otp = as_int(cnt[part]["ot_phylo_extant_species"], 0)
        rows.append(OrderedDict([
            ("partition", part), ("described", as_int(r["described"], 0)),
            ("atlas_C0", as_int(r["c0_units"], 0)), ("atlas_C1", c1),
            ("atlas_C1_tips_enumerated", "yes" if mkeys is not None else "no"),
            ("ot_tips_extant_species",
             as_int(cnt[part]["ot_tips_extant_species"], 0)),
            ("ot_phylo_extant_species", otp),
            ("atlas_C1_over_ot_phylo", round(c1 / otp, 4) if otp else ""),
            ("keys_atlas", len(akeys)), ("keys_ot_phylo", len(o)),
            ("keys_both", len(akeys & o)), ("keys_atlas_only", len(akeys - o)),
            ("keys_ot_only", len(o - akeys)), ("keys_union", len(akeys | o)),
            ("keys_atlas_molecular", len(mkeys) if mkeys is not None else ""),
            ("keys_molecular_and_ot", len(mkeys & o) if mkeys is not None
             else ""),
            ("keys_molecular_union_ot", len(mkeys | o) if mkeys is not None
             else ""),
            ("accepted_molecular", len(macc) if macc is not None else ""),
            ("accepted_molecular_and_ot", len(macc & oa) if macc is not None
             else ""),
            ("accepted_molecular_not_in_ot", len(macc - oa)
             if macc is not None else ""),
        ]))
    den = sum(as_int(r["described"], 0) for r in basis if is_euk(r))
    tot_c1 = round(sum(x["atlas_C1"] for x in rows), NUMERATOR_DP)
    tot_ot = sum(x["ot_phylo_extant_species"] for x in rows)
    euk_ot = as_int(cnt["__Eukaryota"]["ot_phylo_extant_species"], 0)
    sums = OrderedDict((k, sum(x[k] for x in rows if x[k] != ""))
                       for k in ("keys_both", "keys_union", "keys_atlas",
                                 "keys_ot_phylo"))
    enum = [x for x in rows if x["atlas_C1_tips_enumerated"] == "yes"]
    block = OrderedDict([
        ("definition", opentree_comparison.__doc__.split("\n\n")[1].strip()),
        ("release", cnt["__all"]["release"]),
        ("ot_tips_all_life", as_int(cnt["__all"]["ot_tips"], 0)),
        ("ot_phylo_all_life", as_int(cnt["__all"]["ot_phylo"], 0)),
        ("ot_phylo_extant_species_eukaryota", euk_ot),
        ("ot_phylo_extant_species_eukaryota_pct", pct(euk_ot, den)),
        ("n_partitions_compared", len(rows)),
        ("ot_phylo_extant_species_in_atlas_clades", tot_ot),
        ("ot_phylo_extant_species_in_atlas_clades_pct", pct(tot_ot, den)),
        ("atlas_C1_in_compared_partitions", tot_c1),
        ("atlas_C1_in_compared_partitions_pct", pct(tot_c1, den)),
        ("atlas_C1_over_ot_phylo", round(tot_c1 / tot_ot, 4)),
        ("denominator", den),
        ("keys", sums),
        ("molecular_enumerated", OrderedDict([
            ("partitions", [x["partition"] for x in enum]),
            ("keys_atlas_molecular", sum(x["keys_atlas_molecular"]
                                         for x in enum)),
            ("keys_ot_phylo", sum(x["keys_ot_phylo"] for x in enum)),
            ("keys_molecular_and_ot", sum(x["keys_molecular_and_ot"]
                                          for x in enum)),
            ("keys_molecular_union_ot", sum(x["keys_molecular_union_ot"]
                                            for x in enum)),
            ("accepted_rule",
             "R1-14: molecular tips and Open Tree phylogeny-placed species "
             "both resolved to COL26.9 accepted names by the S05 rule "
             "(scripts/match_col.py), compared by accepted-name id, so "
             "synonymy no longer separates the two sets"),
            ("accepted_molecular", sum(x["accepted_molecular"] for x in enum)),
            ("accepted_molecular_and_ot", sum(x["accepted_molecular_and_ot"]
                                              for x in enum)),
            ("accepted_molecular_not_in_ot",
             sum(x["accepted_molecular_not_in_ot"] for x in enum)),
        ])),
        ("per_partition", OrderedDict((x["partition"], OrderedDict(
            (k, x[k]) for k in OT_HEADER if k != "partition")) for x in rows)),
    ])
    return block


def build_opentree_csv(head, warnings):
    rows = [OrderedDict([("partition", p)] + list(v.items()))
            for p, v in head["opentree_P01"]["per_partition"].items()]
    return render_csv(OT_HEADER, rows)


WORKED_FITS = os.path.join(ROOT, "Paper", "inputs", "worked_example_bd.csv")
WORKED_SEL = os.path.join(ROOT, "Paper", "inputs",
                          "worked_example_bd_selection.csv")


WORKED_POLY = os.path.join(ROOT, "Paper", "inputs",
                           "worked_example_bd_polytomy.csv")
WORKED_FORDER = os.path.join(ROOT, "Paper", "inputs",
                             "worked_example_bd_family_order.csv")
WORKED_SETS = os.path.join(ROOT, "Paper", "inputs",
                           "worked_example_bd_sets.csv")
WORKED_SETS_REGISTER = os.path.join(ROOT, "Data",
                                    "worked_example_tree_sets.csv")
# round F2 (F2-04, F2-05, F2-35): the tree sets the fits are repeated on, in
# register order, with the name each takes in macros and in Table S-sets
WORKED_SET_TAG = OrderedDict([
    ("posterior", ("Posterior", "posterior set")),
    ("date_random", ("DateRandom", "date assignment, random dates")),
    ("date_mean", ("DateMean", "date assignment, mean dates")),
    ("prior", ("Prior", "root ages drawn from the prior")),
    ("molecular_only", ("MolOnly", "source's tree of sequenced species")),
])
WORKED_FAMILY_MIN_TIPS = 20   # D23: families this large enter the dispersion
# A15.  The rho identity: two fits on one tree are the same fit when their
# log-likelihoods agree to the 4 dp the script writes and lambda scales by the
# inverse rho ratio to the 6 significant figures it writes.
WORKED_IDENTITY_LL_TOL = 1e-3
WORKED_IDENTITY_RATIO_RTOL = 1e-4
# D27.  The misapplied fit constrains the matching one by lambda*rho >= r*rho'
# (mu >= 0), an inequality on one parameter: the likelihood-ratio statistic is
# a 50:50 mixture of chi2(0) and chi2(1) (Self & Liang 1987), whose 5% critical
# value is qchisq(0.90, 1).
WORKED_POSTERIOR_SET = "posterior"   # Data/worked_example_tree_sets.csv set
WORKED_LR_CUTOFF = 2.705543   # qchisq(0.90, 1): the 5% critical value of one test
WORKED_TEST_LEVEL = 0.05      # the level WORKED_LR_CUTOFF is the critical value for


def worked_lr_cutoff(n_tests):
    """Round F4 (F4-01, D50).  The family of tests is every eligible group, a
    number fixed by the eligibility rule before any fit is read; each test is
    made at the Bonferroni level WORKED_TEST_LEVEL / n_tests, whose critical
    value under the 50:50 chi2(0):chi2(1) mixture is qchisq(1 - 2a, 1) =
    z_{1-a}^2."""
    import statistics as _st
    return _st.NormalDist().inv_cdf(1.0 - WORKED_TEST_LEVEL / n_tests) ** 2
WORKED_EPS_BOUNDARY = 1e-6    # = EPS_BOUNDARY in worked_example_bd.R


def worked_eligible(basis, mol, shipped):
    """The partitions the worked example fits (D23), with their molecular tip
    ids: eukaryotic canonical chronograms in absolute time whose molecular tips
    are enumerated and fewer than all tips.  All are reported; none is
    selected."""
    out = OrderedDict()
    for r in sorted(basis, key=lambda x: x["partition"]):
        p = r["partition"]
        ev = mol.get(p)
        if not (is_euk(r) and r["is_chronogram"] == "yes"
                and r["time_scale"] == "absolute" and ev is not None
                and ev.tip_names is not None and ev.source_file):
            continue
        idmap = shipped.ids.get(p, {})
        ids = {idmap[nm] for nm in ev.tip_names if nm in idmap}
        if len(ids) < as_int(r["labels_unique"], 0):
            out[p] = ids
    return out


def family_sampling(part, mol_ids, shipped):
    """Per CoL family of the canonical tree (families with at least
    WORKED_FAMILY_MIN_TIPS tips), the fraction of its tips that are molecular.
    Uniform sampling of the molecular-pruned tree would make these fractions
    equal; their spread measures how far that assumption fails."""
    fam = defaultdict(lambda: [0, 0])
    for sid, r in shipped.col_rows.get(part, {}).items():
        f = (r.get("col_family") or "").strip()
        if not f:
            continue
        fam[f][0] += 1
        fam[f][1] += 1 if sid in mol_ids else 0
    fr = sorted(m / float(n) for n, m in fam.values()
                if n >= WORKED_FAMILY_MIN_TIPS)
    if not fr:
        return OrderedDict([("n_families", 0)])

    def q(x):
        k = (len(fr) - 1) * x
        lo = int(k)
        hi = min(lo + 1, len(fr) - 1)
        return round(fr[lo] + (fr[hi] - fr[lo]) * (k - lo), FRAC_DP)
    return OrderedDict([
        ("min_tips_per_family", WORKED_FAMILY_MIN_TIPS),
        ("n_families", len(fr)),
        ("q10", q(0.1)), ("median", q(0.5)), ("q90", q(0.9)),
        ("n_families_no_molecular_tip", sum(1 for x in fr if x == 0.0)),
        ("n_families_all_molecular", sum(1 for x in fr if x == 1.0)),
    ])


def worked_identity(rho_t, rho_c, lam_t, lam_c, ll_t, ll_c, eps_c, cutoff):
    """A15, D27.  Classify one same-tree pair of fits (rho = molecular tips,
    rho = full-tree tips)."""
    pred = rho_t / rho_c
    ratio = lam_c / lam_t
    dll = ll_t - ll_c
    ident = (abs(dll) <= WORKED_IDENTITY_LL_TOL
             and abs(ratio - pred) <= WORKED_IDENTITY_RATIO_RTOL * pred)
    lr = 2.0 * dll
    return OrderedDict([
        ("identity_lambda_ratio_predicted", round(pred, 6)),
        ("loglik_rho_tips", ll_t), ("loglik_rho_full", ll_c),
        ("lr_statistic", round(lr, 4)),
        ("classification", "identity" if ident else "boundary"),
        ("epsilon_rho_full_on_boundary", eps_c < WORKED_EPS_BOUNDARY),
        ("rho_full_rejected", (not ident) and lr > cutoff),
    ])


def worked_example(basis, mol, shipped, crown, warnings):
    """S13 (D23).  The worked example's results, read from the outputs of
    Paper/code/analysis/worked_example_bd.R (R, run separately).  The rebuild
    refuses to report them when they were computed on another deposit: the
    eligible set, the tip counts, the crown ages and the time-scale factors
    are re-derived here and must match the script's selection file."""
    if not (os.path.exists(WORKED_FITS) and os.path.exists(WORKED_SEL)):
        die("the S13 outputs are missing; run Rscript "
            "Paper/code/analysis/worked_example_bd.R")
    elig = worked_eligible(basis, mol, shipped)
    by = {r["partition"]: r for r in basis}
    ages = crown["per_partition"]
    sel = read_csv(WORKED_SEL)
    if [x["partition"] for x in sel] != list(elig):
        die("%s was computed on another deposit (eligible %s, now %s); "
            "re-run Paper/code/analysis/worked_example_bd.R"
            % (rel(WORKED_SEL), [x["partition"] for x in sel], list(elig)))
    for x in sel:
        p = x["partition"]
        want = (as_int(by[p]["labels_unique"]), len(elig[p]),
                as_int(by[p]["c0_units"]), as_int(by[p]["described"]),
                as_float(ages[p]["crown_age_ma"]),
                as_float(ages[p]["geometry_root_ma"]))
        got = (as_int(x["labels_unique"]), as_int(x["n_molecular"]),
               as_int(x["c0_units"]), as_int(x["described"]),
               as_float(x["crown_age_ma"]), as_float(x["geometry_root_ma"]))
        if want != got:
            die("%s: %s was computed on another deposit (tips, molecular "
                "tips, C0, described, crown age, geometry root: %s, now %s); re-run "
                "Paper/code/analysis/worked_example_bd.R"
                % (p, rel(WORKED_SEL), got, want))
    fits = read_csv(WORKED_FITS)
    per = OrderedDict()
    for f in fits:
        key = "%s|%s|%s" % (f["partition"], f["tree"], f["rho_convention"])
        per[key] = OrderedDict((k, f[k]) for k in (
            "n_tips", "crown_age_ma", "time_scale_factor", "rho", "lambda",
            "mu", "r", "r_lo", "r_hi", "epsilon", "epsilon_lo", "epsilon_hi",
            "epsilon_ci", "loglik", "n_tips_family_fraction",
            "family_uniform_check", "n_events", "n_polytomy_nodes",
            "n_polytomy_events", "max_children"))

    def g(part, tree, conv, k):
        return as_float(per["%s|%s|%s" % (part, tree, conv)][k])

    def overlap(part, a, b, v):
        return (g(part, a[0], a[1], v + "_lo") <= g(part, b[0], b[1], v + "_hi")
                and g(part, b[0], b[1], v + "_lo")
                <= g(part, a[0], a[1], v + "_hi"))

    rho_effect, imputation, families, family_fit = (
        OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict())
    M, MC, F = ("molecular", "tips"), ("molecular", "full_misapplied"), \
        ("full", "tips")
    for part in elig:
        rho_effect[part] = OrderedDict([
            ("rho_tips", g(part, *M, k="rho")),
            ("rho_full", g(part, *MC, k="rho")),
            ("epsilon_rho_tips", g(part, *M, k="epsilon")),
            ("epsilon_rho_full", g(part, *MC, k="epsilon")),
            ("lambda_ratio_full_over_tips",
             round(g(part, *MC, k="lambda") / g(part, *M, k="lambda"), 4)),
            ("mu_rho_tips", g(part, *M, k="mu")),
            ("mu_rho_full", g(part, *MC, k="mu")),
            ("r_ratio_full_over_tips",
             round(g(part, *MC, k="r") / g(part, *M, k="r"), 6)),
            ("epsilon_intervals_overlap", overlap(part, M, MC, "epsilon")),
        ])
        rho_effect[part].update(worked_identity(
            g(part, *M, k="rho"), g(part, *MC, k="rho"),
            g(part, *M, k="lambda"), g(part, *MC, k="lambda"),
            g(part, *M, k="loglik"), g(part, *MC, k="loglik"),
            g(part, *MC, k="epsilon"), worked_lr_cutoff(len(elig))))
        FM = ("molecular", "family")
        family_fit[part] = OrderedDict(
            [(k, g(part, *FM, k=k)) for k in (
                "lambda", "mu", "r", "r_lo", "r_hi", "epsilon", "epsilon_lo",
                "epsilon_hi", "loglik", "n_tips", "n_tips_family_fraction",
                "family_uniform_check")]
            + [("epsilon_ci", per["%s|molecular|family" % part]["epsilon_ci"]),
               ("epsilon_uniform_tips", g(part, *M, k="epsilon")),
               ("r_uniform_tips", g(part, *M, k="r")),
               ("epsilon_intervals_overlap_uniform",
                overlap(part, M, FM, "epsilon"))])
        rf, rm = g(part, *F, k="r"), g(part, *M, k="r")
        imputation[part] = OrderedDict([
            ("r_full", rf), ("r_molecular", rm),
            ("r_ratio_full_over_molecular", round(rf / rm, 4)),
            ("r_intervals_overlap", overlap(part, F, M, "r")),
            ("epsilon_full", g(part, *F, k="epsilon")),
            ("epsilon_molecular", g(part, *M, k="epsilon")),
            ("epsilon_intervals_overlap", overlap(part, F, M, "epsilon")),
        ])
        families[part] = family_sampling(part, elig[part], shipped)
    polytomy = worked_polytomy(list(elig), per, worked_lr_cutoff(len(elig)))
    tree_sets = worked_sets(list(elig), rho_effect)
    family_order = worked_family_order(list(elig))
    ratios = [v["r_ratio_full_over_molecular"] for v in imputation.values()]
    return OrderedDict([
        ("definition",
         "constant-rate birth-death with uniform sampling fraction rho "
         "(Stadler 2009) on each eligible canonical chronogram rescaled to "
         "its reported crown age (D20): the full tree and the tree pruned to "
         "its molecular tips at rho = tips/described, and the pruned tree at "
         "rho = full-tree tips/described (the full-tree fraction misapplied to the molecular "
         "tree); 95% profile-likelihood intervals, conditional on the single "
         "tree fitted (no topological or dating uncertainty)"),
        ("eligibility",
         "eukaryotic canonical chronogram in absolute time whose molecular "
         "tips are enumerated and fewer than all tips; every eligible "
         "partition is reported, none is selected"),
        ("eligible", list(elig)),
        ("test_family", OrderedDict([
            ("rule", "one likelihood-ratio test per eligible group, on its "
                     "summary tree; the family is every eligible group, fixed "
                     "by the eligibility rule before any fit is read; "
                     "Bonferroni level = 0.05 / n_tests (D50)"),
            ("n_tests", len(elig)),
            ("level", WORKED_TEST_LEVEL / len(elig)),
            ("critical_value", round(worked_lr_cutoff(len(elig)), 6)),
        ])),
        ("identity",
         "the likelihood depends on lambda and rho only through lambda*rho "
         "and r (Stadler 2009), so a fit at rho' is the fit at rho with "
         "lambda scaled by rho/rho' and r and the log-likelihood unchanged, "
         "unless that forces mu < 0; then mu = 0 and the likelihood drops "
         "(the rho-epsilon confounding: Cusimano & Renner 2010 Syst Biol "
         "59:458; Hoehna et al. 2011 Mol Biol Evol 28:2577). Each pair is "
         "classified 'identity' or 'boundary'; at the boundary lr_statistic "
         "= 2 (logL(rho_tips) - logL(rho_full)) and rho_full is rejected when it "
         "exceeds the critical value of the 50:50 chi2(0):chi2(1) mixture at "
         "the Bonferroni level 0.05 / (number of eligible groups) (test_family, "
         "F4-01, D50): "
         "in the parameters of the matching fit, the fit at rho_full is the "
         "same model restricted to epsilon >= 1 - rho_tips/rho_full, an edge "
         "inside (0, 1), so the null distribution at that edge is the mixture "
         "(Chernoff 1954; Self & Liang 1987; D27, D30). Every interval uses "
         "the cutoff qchisq(0.95, 1)/2; where the epsilon estimate is 0 the "
         "interval is [0, upper] (D30)"),
        ("rho_effect_same_tree", rho_effect),
        ("family_sampling_fit", OrderedDict([
            ("definition",
             "the molecular tree with each tip's sampling fraction = accepted "
             "names among the family's molecular tips / the family's accepted "
             "species in COL26.9 (capped at 1; tips with no resolved family "
             "take rho = tips/described), clade-specific likelihood of "
             "FitzJohn, Maddison & Otto 2009; the rho = tips/described fits "
             "are conditional on uniform sampling, which the per-family "
             "fractions contradict. The family and uniform likelihoods carry "
             "different probabilities of sampling the same tips, terms free "
             "of the rates, so they are not compared; the estimates are "
             "reported side by side (D30)"),
            ("per_partition", family_fit)])),
        ("full_vs_molecular", OrderedDict([
            ("note",
             "two different trees, not a rung effect: the difference is what "
             "taxonomically placed tips do to the estimate (Rabosky 2015 "
             "Evolution 69:3207; Chang, Rabosky & Alfaro 2020 Syst Biol "
             "69:602)"),
            ("n_partitions", len(ratios)),
            ("n_r_full_greater", sum(1 for x in ratios if x > 1)),
            ("n_r_full_smaller", sum(1 for x in ratios if x < 1)),
            ("n_r_intervals_overlap",
             sum(1 for v in imputation.values() if v["r_intervals_overlap"])),
            ("per_partition", imputation),
        ])),
        ("molecular_fraction_by_family", families),
        ("polytomies", polytomy),
        ("tree_sets", tree_sets),
        ("family_child_order", family_order),
        ("fits", per),
    ])


def worked_polytomy(elig, per, cutoff):
    """Round F1 (F1-14).  Per eligible group: the branching events of each
    fitted tree that sit at polytomies, and the uniform-sampling fits repeated
    with the polytomies resolved at random (worked_example_bd.R, one row per
    seed): r on each tree, the ratio of r on the full to the pruned tree, the
    extinction fraction of the two pruned fits, and the likelihood-ratio
    statistic of the full tree's fraction on the pruned tree."""
    if not os.path.exists(WORKED_POLY):
        die("%s is missing; run Rscript Paper/code/analysis/worked_example_bd.R"
            % rel(WORKED_POLY))
    rows = read_csv(WORKED_POLY)
    out = OrderedDict()
    for p in elig:
        fit = {}
        for r in rows:
            if r["partition"] == p:
                fit.setdefault(as_int(r["seed"]), {})[
                    (r["tree"], r["rho_convention"])] = r
        if not fit:
            die("%s has no rows for %s" % (rel(WORKED_POLY), p))
        seeds = OrderedDict()
        for sd in sorted(fit):
            f = fit[sd]
            rf = as_float(f[("full", "tips")]["r"])
            rm = as_float(f[("molecular", "tips")]["r"])
            lr = 2.0 * (as_float(f[("molecular", "tips")]["loglik"])
                        - as_float(f[("molecular", "full_misapplied")]["loglik"]))
            em = as_float(f[("molecular", "full_misapplied")]["epsilon"])
            seeds[sd] = OrderedDict([
                ("r_full", rf), ("r_molecular", rm),
                ("r_ratio_full_over_molecular", round(rf / rm, 4)),
                ("epsilon_rho_tips",
                 as_float(f[("molecular", "tips")]["epsilon"])),
                ("epsilon_rho_full", em),
                ("lr_statistic", round(lr, 4)),
                ("boundary", em < WORKED_EPS_BOUNDARY
                 and lr > 2 * WORKED_IDENTITY_LL_TOL),
                ("rejected", em < WORKED_EPS_BOUNDARY
                 and lr > cutoff),
            ])
        ev = OrderedDict()
        for tree in ("full", "molecular"):
            x = per["%s|%s|tips" % (p, tree)]
            ev[tree] = OrderedDict([
                ("n_events", as_int(x["n_events"])),
                ("n_polytomy_events", as_int(x["n_polytomy_events"])),
                ("n_polytomy_nodes", as_int(x["n_polytomy_nodes"])),
                ("max_children", as_int(x["max_children"])),
                ("polytomy_events_pct", pct(as_int(x["n_polytomy_events"]),
                                            as_int(x["n_events"]))),
            ])
        out[p] = OrderedDict([("events", ev), ("resolved", seeds)])
    return out


def worked_sets(elig, eff):
    """Round F2 (F2-04, F2-05, F2-35).  The fits of worked_example_sets.R on
    the tree sets the sources archive (Data/worked_example_tree_sets.csv),
    summarised per set: trees fitted, tips, the extinction fraction of the
    matching fit, the distribution of the likelihood-ratio statistic of the
    full tree's fraction and the trees in which the restriction binds, and,
    where the full tree is fitted too, the ratio of r on the full to r on the
    pruned tree.  The trees are draws from one archived set, not independent
    replicates, so no count of trees above a critical value is formed (F4-03):
    whether the restriction binds, and the statistic, are properties of each
    tree."""
    for f_ in (WORKED_SETS, WORKED_SETS_REGISTER):
        if not os.path.exists(f_):
            die("%s is missing; run Rscript Paper/code/analysis/"
                "worked_example_sets.R" % rel(f_))
    reg = read_csv(WORKED_SETS_REGISTER)
    rows = read_csv(WORKED_SETS)
    out = OrderedDict()
    for g in reg:
        p, st = g["partition"], g["set"]
        if p not in elig:
            die("%s: %s is not an eligible group" % (rel(WORKED_SETS_REGISTER), p))
        if st not in WORKED_SET_TAG:
            die("%s: unknown set %r" % (rel(WORKED_SETS_REGISTER), st))
        rs = [r for r in rows if r["partition"] == p and r["set"] == st]
        if not rs:
            die("%s has no rows for %s %s" % (rel(WORKED_SETS), p, st))
        eps = [as_float(r["epsilon_matching"]) for r in rs]
        lrs = [as_float(r["lr_statistic"]) for r in rs]
        rr = [as_float(r["r_full_tree"]) / as_float(r["r_matching"])
              for r in rs if r["r_full_tree"] not in ("", None)]
        bound = [as_float(r["epsilon_full_rho"]) < WORKED_EPS_BOUNDARY
                 and as_float(r["lr_statistic"]) > 2 * WORKED_IDENTITY_LL_TOL
                 for r in rs]
        out["%s|%s" % (p, st)] = OrderedDict([
            ("partition", p), ("set", st), ("tips", g["tips"]),
            ("n_trees", len(rs)),
            ("n_tips_pruned", sorted({as_int(r["n_tips"]) for r in rs})),
            ("n_tips_full", sorted({as_int(r["n_tips_full_tree"]) for r in rs
                                    if r["n_tips_full_tree"] not in ("", None)})),
            ("max_polytomy_events", max(as_int(r["n_polytomy_events"]) for r in rs)),
            ("epsilon_matching_min", min(eps)),
            ("epsilon_matching_median", median(eps)),
            ("epsilon_matching_max", max(eps)),
            ("lr_min", min(lrs)), ("lr_median", median(lrs)), ("lr_max", max(lrs)),
            ("n_boundary", sum(bound)),
            ("r_ratio_min", min(rr) if rr else None),
            ("r_ratio_max", max(rr) if rr else None),
            ("summary_tree_classification", eff[p]["classification"]),
            ("summary_tree_lr", eff[p]["lr_statistic"]),
            ("summary_tree_epsilon_matching", eff[p]["epsilon_rho_tips"]),
        ])
    return out


def worked_family_order(elig):
    """Round F1 (F1-16).  The family fit with u at each node taken from the
    first and from the last child (worked_example_bd.R)."""
    if not os.path.exists(WORKED_FORDER):
        die("%s is missing; run Rscript Paper/code/analysis/worked_example_bd.R"
            % rel(WORKED_FORDER))
    out = OrderedDict()
    rows = read_csv(WORKED_FORDER)
    for p in elig:
        by = {r["child"]: r for r in rows if r["partition"] == p}
        if set(by) != {"first", "last"}:
            die("%s: %s needs a first and a last row" % (rel(WORKED_FORDER), p))
        a, b = by["first"], by["last"]
        out[p] = OrderedDict([
            ("r_first", as_float(a["r"])), ("r_last", as_float(b["r"])),
            ("epsilon_first", as_float(a["epsilon"])),
            ("epsilon_last", as_float(b["epsilon"])),
            ("loglik_first", as_float(a["loglik"])),
            ("loglik_last_at_first_mle", as_float(a["loglik_last_at_first_mle"])),
            ("r_rel_diff", abs(as_float(b["r"]) - as_float(a["r"]))
             / as_float(a["r"])),
            ("epsilon_abs_diff", abs(as_float(b["epsilon"])
                                     - as_float(a["epsilon"]))),
            ("loglik_abs_diff_at_mle", abs(as_float(a["loglik_last_at_first_mle"])
                                           - as_float(a["loglik"]))),
        ])
    return out


# ==========================================================================
# SUPPLEMENT -- the generated tables and macros of
# Paper/shared/supplementary_body_v2.tex.  Every number the Supplement prints
# is a cell of one of these tables or a macro of supp_macros.tex; each table
# is written twice, as a CSV that verify_atlas.py (SUPPLEMENT) re-derives
# cell by cell, and as the LaTeX body rows rendered from the same rows.
# ==========================================================================

REFERENCE_TREE_SOURCES = os.path.join(ROOT, "Data", "reference_tree_sources.csv")
CANONICAL_SOURCES = os.path.join(ROOT, "Data", "canonical_sources.csv")
CANONICAL_SOURCE_PREP = os.path.join(ROOT, "Data", "canonical_source_prep.csv")
COL_DENOMINATORS = os.path.join(ROOT, "Data", "col_denominators.csv")
NAME_OVERRIDES = os.path.join(STD_DIR, "name_overrides.csv")
NAME_REJECTED = os.path.join(STD_DIR, "name_rejected_targets.csv")
SUPRA_GENERA = os.path.join(STD_DIR, "col269_suprageneric_genera.txt")
GBIF_UNMATCHED = os.path.join(ROOT, "Data", "external", "gbif_unmatched.csv")
WORKED_FITS = os.path.join(ROOT, "Paper", "inputs", "worked_example_bd.csv")
WORKED_SEL = os.path.join(ROOT, "Paper", "inputs",
                          "worked_example_bd_selection.csv")
WORKED_CHECK = os.path.join(ROOT, "Paper", "inputs",
                            "worked_example_bd_crosscheck.csv")
WORKED_FAMILY = os.path.join(ROOT, "Paper", "inputs",
                             "worked_example_bd_family_fractions.csv")
RESTORED_RECORD = os.path.join(ROOT, "Data", "_archive", "2026-10-revision",
                               "corrupted_branch_lengths", "RESTORED.csv")
TREE_REFERENCES = os.path.join(STD_DIR, "tree_references.csv")
VERIFY_SCRIPT = os.path.join(SCRIPT_DIR, "verify_atlas.py")
CROWN_SCRIPT = os.path.join(SCRIPT_DIR, "crown_intervals.py")
SUPP_DIR = os.path.join(ROOT, "Paper", "inputs")

# selection codes of Data/reference_tree_sources.csv, as printed (Table S2)
SUPP_SELECTION_CODE = OrderedDict([
    ("as_archived", "A"), ("release", "R"), ("study_chronogram", "C"),
    ("more_branch_lengths", "B"), ("first_of_set", "F"), ("mcc_of_set", "M"),
])
# kinds of crown interval (Data/crown_intervals.csv), as printed (Table S7)
SUPP_KIND_CODE = OrderedDict([
    ("posterior_sample", "post"), ("bootstrap_dating_sample", "boot"),
    ("node_hpd_annotation", "hpd"), ("node_label_interval", "lab"),
    ("calibration_prior_sample", "prior"), ("taxon_addition_sample", "tax"),
    ("date_assignment_sample", "date"),
])
SUPP_CAL_CODE = OrderedDict([
    ("two_sided", "T"), ("one_sided", "O"), ("none", "N"),
    ("not_checked", "?"),
])
SUPP_MOLSOURCE_CODE = OrderedDict([
    ("source_tree_intersection", "T"), ("source_alignment_intersection", "L"),
    ("per_tip_evidence", "E"), ("asserted_no_imputation", "S"),
])
SUPP_SENS_ROWS = OrderedDict([
    ("reported", "Reported"),
    ("insects_minus", "Insect count $-$\\SuppElasticityPct\\%"),
    ("insects_plus", "Insect count $+$\\SuppElasticityPct\\%"),
    ("col_filled", "Catalogue of Life count for every group"),
    ("col_base", "Catalogue of Life count where the base release has one"),
    ("rep", "Groups with a tree only"),
    ("no_seed", "Seed plants removed (tree and described species)"),
    ("drop_automated_synthesis", "Automated synthesis removed"),
    ("drop_compiled", "Compiled trees removed"),
    ("drop_primary_analysis", "Primary analyses removed"),
    ("retain_only_automated_synthesis", "Automated synthesis only"),
    ("retain_only_compiled", "Compiled trees only"),
    ("retain_only_primary_analysis", "Primary analyses only"),
    ("bird_band_low", "Bird molecular fraction at its lower end"),
    ("union", "Named species of all trees"),
    ("rule_p_counted", "Intervals excluded by rule (p) counted"),
    ("crown_calibrated", "Intervals on a calibrated crown or below a "
                         "two-sided calibration not counted"),
])


def _supp_pct(x):
    """Supplement percentages: the prose rule of macro_pct (D9)."""
    return macro_pct(x)


def _leaf_count_text(text):
    """Leaves of a Newick string: labels that follow '(' or ','."""
    text = re.sub(r"\[[^\]]*\]", "", text)
    return len(re.findall(r"[(,]\s*(?:'[^']*'|[^(),:;\s]+)", text))


def _supp_reference_rows(head, basis):
    reg = {r["filename"]: r for r in read_csv(DATASET_REGISTRY)}
    src = {r["partition"]: r for r in read_csv(REFERENCE_TREE_SOURCES)}
    prep = {r["filename"]: r for r in read_csv(CANONICAL_SOURCE_PREP)}
    rows = []
    for b in basis:
        if not is_euk(b) or as_int(b.get("c0_units"), 0) <= 0:
            continue
        p, fn = b["partition"], b["canonical_tree"]
        s = src.get(p)
        if s is None or s["filename"] != fn:
            die("%s: no row for %s (%s)" % (rel(REFERENCE_TREE_SOURCES), p, fn))
        if s["selection"] not in SUPP_SELECTION_CODE:
            die("%s: %s has selection %r" % (rel(REFERENCE_TREE_SOURCES), p,
                                             s["selection"]))
        if fn in prep and prep[fn]["raw_source"] != s["source_file"]:
            die("%s: %s source %s differs from %s"
                % (rel(REFERENCE_TREE_SOURCES), p, s["source_file"],
                   rel(CANONICAL_SOURCE_PREP)))
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("filename", fn), ("bibkey", reg[fn]["bibkey"]),
            ("source_file", os.path.basename(s["source_file"])),
            ("tree_index", s["tree_index"]),
            ("selection", s["selection"]), ("tree_used", s["tree_used"]),
            ("archive", s["archive"]),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _supp_pruned_rows(head, basis, shipped):
    """Per reference tree: tips of the archived source, tips outside the group
    removed (D7), tips merged into one per species (D4), tips kept, and the
    kept tips whose label names no species (D13)."""
    cs = {r["filename"]: r for r in read_csv(CANONICAL_SOURCES)}
    removed = Counter(r["filename"] for r in read_csv(PRUNED_TIPS))
    rows = []
    for b in basis:
        if not is_euk(b) or as_int(b.get("c0_units"), 0) <= 0:
            continue
        p, fn = b["partition"], b["canonical_tree"]
        with open(os.path.join(ROOT, cs[fn]["source"]), encoding="utf-8",
                  errors="replace") as fh:
            n_src = _leaf_count_text(fh.read())
        tips = as_int(b["labels_unique"], 0)
        merged = n_src - removed[fn] - tips
        if merged < 0:
            die("%s: source %d - removed %d < tips %d"
                % (p, n_src, removed[fn], tips))
        open_n = sum(1 for nm, sid in shipped.ids[p].items()
                     if not SPECIES_KEY(nm, sid)[1])
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("source_tips", n_src), ("removed_outside", removed[fn]),
            ("merged_same_species", merged), ("tips", tips),
            ("unresolved_labels_kept", open_n),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _override_counts():
    """Round F1 (F1-06): per group, the resolved labels whose query is the
    binomial of Data/standardized/name_overrides.csv, not the label's own
    species key."""
    n = Counter()
    for r in read_csv(COL_RESOLUTION):
        if r.get("query_source") == "override" \
                and r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
            n[r["partition"]] += 1
    return n


def _supp_names_rows(head, basis):
    nr = head["name_resolution_S05"]["per_partition"]
    ovr = _override_counts()
    gb = head["sensitivities"]["c0_gbif_unmatched"]["per_partition"]
    by = {r["partition"]: r for r in basis}
    rows = []
    for p, v in nr.items():
        if not v.get("eukaryotic"):
            continue
        mt = v["match_types"]
        syn = sum(n for k, n in mt.items() if k.startswith("synonym"))
        acc = sum(n for k, n in mt.items()
                  if k.startswith("accepted") or k.startswith("provisional"))
        unres = v["labels"] - v["resolved_labels"]
        if acc + syn != v["resolved_labels"]:
            die("Supplement S3: %s match types do not sum to its resolved "
                "labels" % p)
        rows.append(OrderedDict([
            ("partition", p), ("category", by[p]["category"]),
            ("labels", v["labels"]), ("accepted", acc), ("synonym", syn),
            ("unresolved", unres),
            ("unresolved_pct", fnum(pct(unres, v["labels"]), 1)),
            ("named_species", v["distinct_accepted_names"]),
            ("merged", v["resolved_labels"] - v["distinct_accepted_names"]),
            ("gbif_added", (gb.get(p) or {}).get("gbif_species_added", 0)),
            ("override", ovr[p]),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _supp_molecular_rows(head, basis):
    ev = head["molecular_evidence"]
    prod = head["sensitivities"]["c1_definition"]["per_partition"]
    molonly = _molecular_only_trees()
    modat = _molecular_only_dating()
    nmol = _molecular_tip_counts()
    rows = []
    for b in basis:
        if not is_euk(b) or as_int(b.get("c0_units"), 0) <= 0:
            continue
        p = b["partition"]
        e = ev[p]
        tips = as_int(b["labels_unique"], 0)
        mt = nmol[p] if b["c1_basis"] == "enumerated" else tips
        srcf = b.get("molecular_source_file") or ""
        basis_ = b["f_basis"]
        if basis_ not in SUPP_MOLSOURCE_CODE:
            die("Supplement S4: %s has f_basis %r" % (p, basis_))
        exact = e.get("exact_label_intersection")
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("f_basis", basis_),
            ("source_file", os.path.basename(srcf) if srcf else ""),
            ("tips", tips), ("molecular_tips", mt),
            ("molecular_fraction", "%.3f" % (mt / tips)),
            ("exact_label_tips", "" if exact is None else int(exact)),
            ("molecular_species", int(round(as_float(b["c1_units"])))),
            ("product", fnum(prod[p]["c1_product"], 1)),
            ("f_low", "%.3f" % as_float(b["f_low"])),
            ("f_high", "%.3f" % as_float(b["f_high"])),
            ("molecular_only_tree", "yes" if p in molonly else "no"),
            ("molecular_only_dated", ("yes" if modat.get(p) == "time" else "no")
             if p in molonly else ""),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _supp_time_rows(head, basis):
    sd = _study_dating()
    rows = []
    for b in basis:
        if not is_euk(b) or as_int(b.get("c0_units"), 0) <= 0:
            continue
        p, fn = b["partition"], b["canonical_tree"]
        h = ROOT_HEIGHT.get(fn)
        if h is None:
            die("Supplement S5: no measured root height for %s" % fn)
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("roottotip_cv", "%.4f" % as_float(b["roottotip_cv"])),
            ("ultrametric", "yes" if as_float(b["roottotip_cv"])
             <= CHRONOGRAM_TOLERANCE else "no"),
            ("root_height", "%.4g" % h),
            ("time_scale", b["time_scale"] or ""),
            ("study_dating", sd.get(p) or (
                "absolute_deposited" if b["is_chronogram"] == "yes"
                else "relative_only" if b["time_scale"] == "relative" else "")),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _supp_crown_rows(head, basis):
    ca = head["crown_ages"]["per_partition"]
    civ = crown_intervals()
    calreg = calibration_register()
    cat = {r["partition"]: r["category"] for r in basis}
    rows = []
    for p, c in ca.items():
        kind = c["interval_kind"]
        iv = civ.get(p) or {}
        if not kind:
            status = "none"
        elif c["counts_for_c1b"] == "yes":
            status = "counted"
        else:
            status = c["c1b_exclusion"]
        pc, pw = iv.get("prior_containment") or "", iv.get("prior_width_ratio") or ""
        rows.append(OrderedDict([
            ("partition", p), ("category", cat[p]),
            ("crown_clade", c["crown_clade"]),
            ("crown_age_ma", "%.1f" % as_float(c["crown_age_ma"])),
            ("age_statistic", c["age_statistic"]),
            ("interval_lo_ma", "%.1f" % as_float(c["interval_lo_ma"])
             if kind else ""),
            ("interval_hi_ma", "%.1f" % as_float(c["interval_hi_ma"])
             if kind else ""),
            ("interval_width_myr", "%.3g" % as_float(c["interval_width_myr"])
             if kind else ""),
            ("interval_kind", kind), ("n_trees", c["n_trees"] if kind else ""),
            ("calibration", c["calibration"] if kind else ""),
            ("prior_lo_ma", "%.1f" % as_float(iv["prior_lo_ma"])
             if iv.get("prior_lo_ma") else ""),
            ("prior_hi_ma", "%.1f" % as_float(iv["prior_hi_ma"])
             if iv.get("prior_hi_ma") else ""),
            ("prior_containment_pct", "%.0f" % (100 * as_float(pc))
             if pc else ""),
            ("prior_width_ratio", "%.2f" % as_float(pw) if pw else ""),
            # F3-14: where the interval and its centre lie in the specified
            # prior density (Data/calibration_register.csv prior_density)
            ("prior_q_lo", iv.get("prior_q_lo") or "" if kind else ""),
            ("prior_q_center", iv.get("prior_q_center") or "" if kind else ""),
            ("prior_q_hi", iv.get("prior_q_hi") or "" if kind else ""),
            # F3-19: the calibration on the node above the crown
            ("parent_calibration", (calreg.get(p, {}).get(
                "parent_calibration") or "") if kind else ""),
            # F4-09 (D52): the calibrations below the crown, and whether the
            # range rule (p) used is the crown's own prior or the band that
            # the calibrations around an uncalibrated crown allow
            ("child_calibration", (calreg.get(p, {}).get(
                "child_calibration") or "") if kind else ""),
            ("prior_range", iv.get("prior_range") or "" if kind else ""),
            ("status", status),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


def _rung_pcts(rungs, den=None):
    out = OrderedDict()
    for k in LEVEL_ORDER:
        e = rungs[k]
        n = e.get("numerator_union", e.get("numerator_sum", e.get("numerator")))
        d = den if den is not None else e["denominator"]
        out[k] = pct(n, d)
    return out


def _supp_sensitivity_rows(head, basis):
    """The alternatives of Supplement S7, each as the four level percentages
    and the denominator they are formed over."""
    rg, s = head["rungs"], head["sensitivities"]
    D = head["denominators"]["eukaryotic"]["described"]
    ds = s["denominator"]
    el = ds["elasticity_demonstration"]["results"]
    f = ds["elasticity_demonstration"]["fraction"]
    la = s["leave_assembly_class_out"]
    dom = head["disclosures"]["dominance"]
    dom_d = sum(as_int(r["described"], 0) for r in basis
                if r["partition"] == dom["partition"])
    birds = s["f_mol_bands"]
    rows = []

    def add(key, den, pcts):
        rows.append(OrderedDict(
            [("alternative", key), ("denominator", int(round(den)))]
            + [(k, "" if pcts.get(k) is None else _supp_pct(pcts[k]))
               for k in LEVEL_ORDER]))

    add("reported", D, _rung_pcts(rg))
    add("insects_minus", el["minus_20_pct"]["denominator"],
        _rung_pcts(el["minus_20_pct"]["rungs"]))
    add("insects_plus", el["plus_20_pct"]["denominator"],
        _rung_pcts(el["plus_20_pct"]["rungs"]))
    col = ds["col_single_authority"]
    add("col_filled", col["filled"]["denominator"],
        _rung_pcts(col["filled"]["rungs"]))
    add("col_base", col["base_only"]["denominator"],
        _rung_pcts(col["base_only"]["rungs"]))
    rep = ds["represented_partitions_only"]
    add("rep", rep["described"], _rung_pcts(rep["rungs"]))
    add("no_seed", D - dom_d, OrderedDict(
        (k, pct(dom["%s_total" % k] - dom["contribution_to_%s" % k],
                D - dom_d)) for k in LEVEL_ORDER))
    for key in ("drop_automated_synthesis", "drop_compiled",
                "drop_primary_analysis", "retain_only_automated_synthesis",
                "retain_only_compiled", "retain_only_primary_analysis"):
        add(key, D, _rung_pcts(la[key]["rungs"], D))
    if len(birds) != 1:
        die("Supplement S7: expected one molecular-fraction band, found %d"
            % len(birds))
    band = next(iter(birds.values()))
    add("bird_band_low", D, _rung_pcts(band["at_f_low"], D))
    u = head["union_bound_S09"]
    add("union", D, OrderedDict([("C0", pct(u["C0_union"], D)),
                                 ("C1", pct(u["C1_union_high"], D)),
                                 ("C1d", None), ("C1b", None)]))
    for key, sk in (("rule_p_counted", "c1b_rule_p_counted"),
                    ("crown_calibrated", "c1b_crown_calibrated")):
        e = s[sk]["rungs"]["C1b"]
        add(key, D, OrderedDict([("C0", None), ("C1", None), ("C1d", None),
                                 ("C1b", pct(e["numerator_union"], D))]))
    return rows


def _supp_denominator_rows(head, basis):
    col = {r["partition"]: r for r in read_csv(COL_DENOMINATORS)}
    uni = head["union_bound_S09"]["per_partition"]
    rows = []
    for b in basis:
        if not is_euk(b):
            continue
        p = b["partition"]
        d = as_int(b["described"], 0)
        c = col.get(p)
        if c is None:
            die("Supplement S7: %s has no row in %s" % (p, rel(COL_DENOMINATORS)))
        n0 = as_int(b.get("c0_units"), 0)
        cf = as_int(c["col_filled"], 0)
        u = uni.get(p)
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]), ("described", d),
            ("described_source", b["described_source"]),
            ("col_count", cf),
            ("col_filled_from", c["filled_classes"] or ""),
            ("named_pct", _supp_pct(pct(n0, d)) if n0 else ""),
            ("named_col_pct", _supp_pct(pct(n0, cf)) if n0 and cf else ""),
            ("named_all_trees_pct", _supp_pct(pct(u["C0_union"], d))
             if u else ""),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             -x["described"], x["partition"]))
    return rows


def _supp_opentree_rows(head, basis):
    cat = {r["partition"]: r["category"] for r in basis}
    rows = []
    for r in csv.DictReader(io.StringIO(build_opentree_csv(head, []))):
        if r["atlas_C1_tips_enumerated"] != "yes":
            continue
        am = as_int(r["accepted_molecular"], 0)
        both = as_int(r["accepted_molecular_and_ot"], 0)
        km = as_int(r["keys_atlas_molecular"], 0)
        kb = as_int(r["keys_molecular_and_ot"], 0)
        rows.append(OrderedDict([
            ("partition", r["partition"]), ("category", cat[r["partition"]]),
            ("ot_phylo_species", as_int(r["ot_phylo_extant_species"], 0)),
            ("molecular_species", am), ("both", both),
            ("molecular_not_ot", as_int(r["accepted_molecular_not_in_ot"], 0)),
            ("both_pct", fnum(pct(both, am), 1) if am else ""),
            ("molecular_keys", km),
            ("ot_keys", as_int(r["keys_ot_phylo"], 0)), ("keys_both", kb),
            ("keys_both_pct", fnum(pct(kb, km), 1) if km else ""),
            ("keys_both_ot_pct", fnum(pct(kb, as_int(r["keys_ot_phylo"], 0)), 1)
             if as_int(r["keys_ot_phylo"], 0) else ""),
        ]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             x["partition"]))
    return rows


SUPP_FIT_ORDER = (("full", "tips"), ("molecular", "tips"),
                  ("molecular", "full_misapplied"), ("molecular", "family"))


def _sig(x, n=3):
    x = float(x)
    if x == 0 or abs(x) < 1e-9:
        return "0"
    return "%.*g" % (n, x)


def _supp_fit_rows(head):
    fits = read_csv(WORKED_FITS)
    by = {(r["partition"], r["tree"], r["rho_convention"]): r for r in fits}
    rows = []
    for p in head["worked_example_S13"]["eligible"]:
        for tree, conv in SUPP_FIT_ORDER:
            r = by.get((p, tree, conv))
            if r is None:
                die("Supplement S9: no fit %s/%s/%s" % (p, tree, conv))
            rows.append(OrderedDict([
                ("partition", p), ("tree", tree), ("rho_convention", conv),
                ("n_tips", as_int(r["n_tips"])),
                ("rho", "" if conv == "family" else rho_fmt(as_float(r["rho"]))),
                ("lambda", _sig(r["lambda"])), ("mu", _sig(r["mu"])),
                ("r", "%.4f" % as_float(r["r"])),
                ("r_lo", "%.4f" % as_float(r["r_lo"])),
                ("r_hi", "%.4f" % as_float(r["r_hi"])),
                ("epsilon", "%.3f" % as_float(r["epsilon"])),
                ("epsilon_lo", "%.3f" % as_float(r["epsilon_lo"])),
                ("epsilon_hi", "%.3f" % as_float(r["epsilon_hi"])),
                ("loglik", "%.2f" % as_float(r["loglik"])),
                ("n_events", as_int(r["n_events"])),
                ("n_polytomy_events", as_int(r["n_polytomy_events"])),
            ]))
    return rows


def _supp_polytomy_rows(head):
    """Round F1 (F1-14): the uniform-sampling estimates as fitted and, over
    the seeds of worked_example_bd.R, with the polytomies resolved at random."""
    w = head["worked_example_S13"]
    po, eff = w["polytomies"], w["rho_effect_same_tree"]
    fvm = w["full_vs_molecular"]["per_partition"]
    rows = []

    def rng(vals, fmt):
        lo, hi = fmt % min(vals), fmt % max(vals)
        return lo if lo == hi else "%s--%s" % (lo, hi)
    for p in w["eligible"]:
        rs = list(po[p]["resolved"].values())
        e = eff[p]
        rows.append(OrderedDict([
            ("partition", p),
            ("events_full", po[p]["events"]["full"]["n_events"]),
            ("polytomy_events_full",
             po[p]["events"]["full"]["n_polytomy_events"]),
            ("events_pruned", po[p]["events"]["molecular"]["n_events"]),
            ("polytomy_events_pruned",
             po[p]["events"]["molecular"]["n_polytomy_events"]),
            ("r_ratio", "%.2f" % fvm[p]["r_ratio_full_over_molecular"]),
            ("r_ratio_resolved", rng([v["r_ratio_full_over_molecular"]
                                      for v in rs], "%.2f")),
            ("epsilon_matching", "%.2f" % e["epsilon_rho_tips"]),
            ("epsilon_matching_resolved", rng([v["epsilon_rho_tips"]
                                               for v in rs], "%.2f")),
            ("lr", "%.2f" % e["lr_statistic"]),
            ("lr_resolved", rng([v["lr_statistic"] for v in rs], "%.2f")),
            ("outcome_resolved", "mu0_rejected" if all(v["rejected"] for v in rs)
             else "mu0" if all(v["boundary"] for v in rs)
             else "rescaled" if not any(v["boundary"] for v in rs)
             else "mixed"),
        ]))
    return rows


def _supp_family_rows(head):
    fam = read_csv(WORKED_FAMILY)
    fits = {(r["partition"], r["tree"], r["rho_convention"]): r
            for r in read_csv(WORKED_FITS)}
    rows = []
    for p in head["worked_example_S13"]["eligible"]:
        fs = sorted(as_float(r["f"]) for r in fam if r["partition"] == p)
        if not fs:
            die("Supplement S9: no family fractions for %s" % p)
        fm = fits[(p, "molecular", "family")]
        ft = fits[(p, "molecular", "tips")]
        n = len(fs)
        med = fs[n // 2] if n % 2 else 0.5 * (fs[n // 2 - 1] + fs[n // 2])
        fr = [r for r in fam if r["partition"] == p and r["n_col_species"]]
        pooled = (sum(as_int(r["n_accepted_names"]) for r in fr)
                  / float(sum(as_int(r["n_col_species"]) for r in fr)))
        fo = head["worked_example_S13"]["family_child_order"][p]
        rows.append(OrderedDict([
            ("partition", p), ("families", n),
            ("tips_with_family", as_int(fm["n_tips_family_fraction"])),
            ("tips", as_int(ft["n_tips"])),
            ("f_min", "%.2f" % fs[0]), ("f_median", "%.2f" % med),
            ("f_max", "%.2f" % fs[-1]),
            ("f_pooled", rho_fmt(pooled)),
            ("rho_uniform", rho_fmt(as_float(ft["rho"]))),
            ("r_first", "%.5f" % fo["r_first"]),
            ("r_last", "%.5f" % fo["r_last"]),
            ("epsilon_first", "%.4f" % fo["epsilon_first"]),
            ("epsilon_last", "%.4f" % fo["epsilon_last"]),
        ]))
        if abs(as_float(fm["family_uniform_check"])) > WORKED_UNIFORM_TOL:
            die("Supplement S9: %s family likelihood at equal fractions differs "
                "from the uniform one by %s" % (p, fm["family_uniform_check"]))
    return rows


def _verify_check_list():
    """(id, title) of every named check of verify_atlas.py, read from its
    CHECKS list without importing it (verify never imports the rebuild, and
    the rebuild does not run verify)."""
    import ast
    with open(VERIFY_SCRIPT, encoding="utf-8") as fh:
        src = fh.read()
    m = re.search(r"^CHECKS = (\[.*?^\])", src, re.M | re.S)
    if not m:
        die("cannot read the CHECKS list of %s" % rel(VERIFY_SCRIPT))
    return ast.literal_eval(m.group(1))


def _crown_rule_parameters():
    """QUANTILES and the D5/D22/D25 thresholds of scripts/crown_intervals.py,
    read from its source (it is not imported: it reads deposits at import)."""
    import ast
    with open(CROWN_SCRIPT, encoding="utf-8") as fh:
        src = fh.read()
    out = {}
    for name in ("QUANTILES", "DEGENERATE_WIDTH_MYR", "BOUND_TOL_MYR",
                 "PRIOR_CONTAINMENT_MIN", "PRIOR_WIDTH_RATIO_MIN"):
        m = re.search(r"^%s = ([^#\n]+)" % name, src, re.M)
        if not m:
            die("cannot read %s from %s" % (name, rel(CROWN_SCRIPT)))
        out[name] = ast.literal_eval(m.group(1).strip())
    return out


def _restored_counts():
    """Shipped trees whose branch lengths were restored from their
    pre-relabelling copies, and the branch lengths and internal labels that
    changed (Data/_archive/2026-10-revision/corrupted_branch_lengths)."""
    shipped = set(fn for fn in os.listdir(STD_TREES) if fn.endswith(".nwk"))
    n_tree = n_bl = n_lab = 0
    for r in read_csv(RESTORED_RECORD):
        f = r["file"]
        if not f.startswith("Data/standardized/trees/") or r["changed"] != "yes":
            continue
        if os.path.basename(f) not in shipped:
            continue
        n_tree += 1
        n_bl += as_int(r["branch_lengths_restored"], 0)
        n_lab += as_int(r["internal_labels_restored"], 0)
    return n_tree, n_bl, n_lab, len(shipped)


def build_supp_macros(head, basis, shipped, warnings):
    """Paper/inputs/supp_macros.tex: the quantities the Supplement's prose
    quotes that the main-text macro files do not hold."""
    m = OrderedDict()
    ref = _supp_reference_rows(head, basis)
    for code, letter in SUPP_SELECTION_CODE.items():
        m["SuppSel" + crown_macro_name(code.replace("_", " "))] = macro_int(
            sum(1 for r in ref if r["selection"] == code))
    for code in SUPP_SELECTION_CODE:
        m["SuppSelNames" + crown_macro_name(code.replace("_", " "))] = \
            name_list(sorted(r["partition"] for r in ref
                             if r["selection"] == code)) or "none"
    pr = _supp_pruned_rows(head, basis, shipped)
    m["SuppSourceTipsN"] = macro_int(sum(r["source_tips"] for r in pr))
    m["SuppRemovedN"] = macro_int(sum(r["removed_outside"] for r in pr))
    m["SuppRemovedTrees"] = macro_int(sum(1 for r in pr
                                          if r["removed_outside"]))
    m["SuppMergedN"] = macro_int(sum(r["merged_same_species"] for r in pr))
    m["SuppMergedTrees"] = macro_int(sum(1 for r in pr
                                         if r["merged_same_species"]))
    m["SuppOpenLabelsN"] = macro_int(sum(r["unresolved_labels_kept"]
                                         for r in pr))
    m["SuppTipsN"] = macro_int(sum(r["tips"] for r in pr))
    m["SuppOverridesN"] = macro_int(len(read_csv(NAME_OVERRIDES)))
    m["SuppRejectedTargetsN"] = macro_int(len(read_csv(NAME_REJECTED)))
    with open(SUPRA_GENERA, encoding="utf-8") as fh:
        m["SuppSupraGeneraN"] = macro_int(sum(
            1 for ln in fh if ln.strip() and not ln.startswith("#")))
    nm = _supp_names_rows(head, basis)
    m["SuppSynonymLabelsN"] = macro_int(sum(r["synonym"] for r in nm))
    m["SuppMergedNamesN"] = macro_int(sum(r["merged"] for r in nm))
    # round F2 (F2-17): register rows whose binomial resolves, the rest
    m["SuppOverridesResolvedN"] = macro_int(sum(r["override"] for r in nm))
    m["SuppOverridesUnresolvedN"] = macro_int(
        len(read_csv(NAME_OVERRIDES)) - sum(r["override"] for r in nm))
    g = head["sensitivities"]["c0_gbif_unmatched"]
    m["SuppGbifQueriedN"] = macro_int(g["unmatched_tips"])
    m["SuppGbifResolvedN"] = macro_int(g["unmatched_tips_resolved_by_gbif"])
    mo = _supp_molecular_rows(head, basis)
    for code in SUPP_MOLSOURCE_CODE:
        m["SuppMolBasis" + crown_macro_name(code.replace("_", " "))] = \
            macro_int(sum(1 for r in mo if r["f_basis"] == code))
    m["SuppExactZeroNames"] = name_list(sorted(
        r["partition"] for r in mo if r["exact_label_tips"] == 0)) or "none"
    band = next(iter(head["sensitivities"]["f_mol_bands"].values()))
    m["SuppBandLow"] = "%.3f" % band["f_low"]
    m["SuppBandHigh"] = "%.3f" % band["f_high"]
    ts = _supp_time_rows(head, basis)
    chron = [as_float(r["roottotip_cv"]) for r in ts if r["ultrametric"] == "yes"]
    other = [as_float(r["roottotip_cv"]) for r in ts if r["ultrametric"] == "no"]
    m["SuppMaxCVUltrametric"] = "%.4f" % max(chron)
    m["SuppMinCVOther"] = "%.4f" % min(other)
    m["SuppNUltrametric"] = macro_int(len(chron))
    cr = _crown_rule_parameters()
    q = cr["QUANTILES"]
    m["SuppCrownIntervalPct"] = "%.0f" % (100 * (q[-1] - q[0]))
    for key, nm_ in (("degenerate", "Degenerate"), ("at_bound", "AtBound"),
                     ("prior_reproduced", "Prior"),
                     ("prior_sampled", "PriorSampled"), ("kind", "Kind")):
        m["SuppCrown" + nm_] = name_list(sorted(
            p for p, c in head["crown_ages"]["per_partition"].items()
            if c["c1b_exclusion"] == key)) or "none"
    m["SuppCrownNone"] = name_list(sorted(
        p for p, c in head["crown_ages"]["per_partition"].items()
        if not c["interval_kind"]))
    ds = head["sensitivities"]["denominator"]
    m["SuppElasticityPct"] = "%.0f" % (100 * ds["elasticity_demonstration"]
                                       ["fraction"])
    m["SuppInsectSharePct"] = macro_pct(ds["largest_single_count"]
                                        ["share_of_denominator_pct"])
    m["SuppColBaseKeepN"] = macro_int(len(ds["col_single_authority"]
                                          ["base_only"]
                                          ["partitions_keeping_primary"]))
    la = head["sensitivities"]["leave_assembly_class_out"]["classes_present"]
    m["SuppAssemblyPrimaryN"] = macro_int(la[ASSEMBLY_PRIMARY])
    m["SuppAssemblyCompiledN"] = macro_int(la[ASSEMBLY_COMPILED])
    m["SuppAssemblyAutomatedN"] = macro_int(la[ASSEMBLY_AUTOMATED])
    lr = head["sensitivities"]["leave_assembly_class_out"]
    m["SuppAssemblyCompiledNames"] = name_list(
        lr["retain_only_" + ASSEMBLY_COMPILED]["retained_partitions"])
    m["SuppAssemblyAutomatedNames"] = name_list(
        lr["retain_only_" + ASSEMBLY_AUTOMATED]["retained_partitions"])
    u = head["union_bound_S09"]
    m["SuppUnionTreesN"] = macro_int(u["n_noncanonical_trees"])
    fr = _supp_family_rows(head)
    m["SuppFamilyMinTips"] = macro_int(min(r["tips_with_family"] for r in fr))
    chk = read_csv(WORKED_CHECK)
    m["SuppCrossCheckN"] = macro_int(len(chk))
    w = head["worked_example_S13"]
    m["SuppPolySeedsN"] = macro_int(len(next(iter(
        w["polytomies"].values()))["resolved"]))
    fo = w["family_child_order"].values()
    m["SuppFamOrderMaxRDiffPct"] = "%.1f" % (100 * max(v["r_rel_diff"] for v in fo))
    m["SuppFamOrderMaxEpsDiff"] = "%.4f" % max(v["epsilon_abs_diff"] for v in fo)
    m["SuppFamOrderMaxDLL"] = "%.1f" % max(v["loglik_abs_diff_at_mle"] for v in fo)
    # each rate's difference relative to the speciation rate, the scale of
    # both rates: a small extinction rate is set by the likelihood no more
    # precisely than that scale (round F1, a refit with phytools 2.5.2)
    rd = 0.0
    for r in chk:
        lam = max(as_float(r["this_lambda"]), as_float(r["phytools_lambda"]))
        for a, b in (("this_lambda", "phytools_lambda"),
                     ("this_mu", "phytools_mu")):
            x, y = as_float(r[a]), as_float(r[b])
            if lam > 1e-6:
                rd = max(rd, abs(x - y) / lam)
    if rd > WORKED_CROSSCHECK_RTOL:
        die("Supplement S9: rates differ from phytools by %.2g (relative)" % rd)
    m["SuppCrossCheckTolPct"] = "%g" % (100 * WORKED_CROSSCHECK_RTOL)
    m["SuppUniformCheckTol"] = r"$10^{%d}$" % round(math.log10(WORKED_UNIFORM_TOL))
    m["SuppProfileCutoff"] = "%.2f" % (_chisq1_quantile(WORKED_INTERVAL_LEVEL)
                                       / 2.0)
    with open(os.path.join(SCRIPT_DIR, "col_denominators.py"),
              encoding="utf-8") as fh:
        hm = re.search(r"^HOLE_FRACTION = ([0-9.]+)", fh.read(), re.M)
    if not hm:
        die("cannot read HOLE_FRACTION from scripts/col_denominators.py")
    m["SuppColHolePct"] = "%.0f" % (100 * float(hm.group(1)))
    with open(os.path.join(ROOT, "Paper", "code", "analysis",
                           "worked_example_bd.R"), encoding="utf-8") as fh:
        cm = re.search(r"^CROSSCHECK_MAX_TIPS <- ([0-9]+)", fh.read(), re.M)
    if not cm:
        die("cannot read CROSSCHECK_MAX_TIPS from worked_example_bd.R")
    m["SuppCrossCheckMaxTips"] = macro_int(int(cm.group(1)))
    n_tree, n_bl, n_lab, n_ship = _restored_counts()
    m["SuppRestoredTreesN"] = macro_int(n_tree)
    m["SuppRestoredBranchLengthsN"] = macro_int(n_bl)
    m["SuppShippedTreesN"] = macro_int(n_ship)
    m["SuppVerifyChecksN"] = macro_int(len(_verify_check_list()))
    vr = Counter((r.get("rule") or "").strip()
                 for r in read_csv(TREE_REFERENCES))
    m["SuppRefBuiltN"] = macro_int(vr.get("built_from_source", 0))
    m["SuppRefRelabelN"] = macro_int(vr.get("tokens_equal_after_tip_relabel", 0))
    m["SuppRefEdgeN"] = macro_int(vr.get("edge_length_multiset_equal", 0))
    lines = ["% generated by scripts/rebuild_atlas.py (build_supp_macros): the",
             "% quantities quoted in Paper/shared/supplementary_body_v2.tex"]
    lines += [r"\newcommand{\%s}{%s}" % (k, v) for k, v in m.items()]
    return "\n".join(lines) + "\n"


# qchisq(1 - alpha, 1) = z_{1 - alpha/2}^2: the profile-likelihood cutoff of
# every worked-example interval is half of it (worked_example_bd.R, D30); the
# level is the rule parameter WORKED_INTERVAL_LEVEL.
WORKED_INTERVAL_LEVEL = 0.95
# Agreement the Supplement states, as tolerances: the differences themselves
# are rounding noise that changes between machines (a refit on another host
# moved the largest cross-check difference from 0.003 % to 0.002 % and the
# family-uniform differences by orders of magnitude below 1e-9).
WORKED_UNIFORM_TOL = 1e-6        # family log-likelihood at equal f vs uniform
WORKED_CROSSCHECK_RTOL = 1e-4    # rates against phytools::fit.bd, relative


def _chisq1_quantile(level):
    from statistics import NormalDist
    return NormalDist().inv_cdf(0.5 + level / 2.0) ** 2


def _supp_tex_rows(rows, render, ncol, category_key="category"):
    out, prev = [], None
    for r in rows:
        if category_key and r[category_key] != prev:
            if prev is not None:
                out.append(r"\midrule")
            out.append(r"\multicolumn{%d}{@{}l}{\textit{%s}} \\"
                       % (ncol, tex_escape(r[category_key])))
            prev = r[category_key]
        cells = render(r)
        if len(cells) != ncol:
            die("supplement table row has %d cells, expected %d"
                % (len(cells), ncol))
        out.append(" & ".join(cells) + r" \\")
    return "\n".join(out) + "\n"


def _dash(x):
    return "---" if x in ("", None) else x


SUPP_TABLES = OrderedDict()


def supp_table(name, header):
    def deco(fn):
        SUPP_TABLES[name] = (header, fn)
        return fn
    return deco


def build_supp_csv(name, head, basis, shipped):
    header, fn = SUPP_TABLES[name]
    return render_csv(header, fn(head, basis, shipped, tex=False))


def build_supp_tex(name, head, basis, shipped):
    header, fn = SUPP_TABLES[name]
    return fn(head, basis, shipped, tex=True)


def _archive_url(a):
    """'doi:10.x/y' -> https://doi.org/10.x/y; a bare host path -> https://."""
    return ("https://doi.org/" + a[len("doi:"):] if a.startswith("doi:")
            else "https://" + a)


@supp_table("supp_reference_trees", ["partition", "category", "filename",
                                     "bibkey", "source_file", "tree_index",
                                     "selection", "tree_used", "archive"])
def _t_reference(head, basis, shipped, tex):
    rows = _supp_reference_rows(head, basis)
    if not tex:
        return rows

    def render(r):
        return [tex_escape(r["partition"]), r"\citet{%s}" % r["bibkey"],
                tex_escape(r["tree_used"]), SUPP_SELECTION_CODE[r["selection"]],
                r"\url{%s}" % _archive_url(r["archive"])]
    return _supp_tex_rows(rows, render, 5)


@supp_table("supp_pruned", ["partition", "category", "source_tips",
                            "removed_outside", "merged_same_species", "tips",
                            "unresolved_labels_kept"])
def _t_pruned(head, basis, shipped, tex):
    rows = _supp_pruned_rows(head, basis, shipped)
    if not tex:
        return rows
    return _supp_tex_rows(rows, lambda r: [
        tex_escape(r["partition"]), thousands(r["source_tips"]),
        thousands(r["removed_outside"]), thousands(r["merged_same_species"]),
        thousands(r["tips"]), thousands(r["unresolved_labels_kept"])], 6)


@supp_table("supp_names", ["partition", "category", "labels", "accepted",
                           "synonym", "unresolved", "unresolved_pct",
                           "named_species", "merged", "gbif_added",
                           "override"])
def _t_names(head, basis, shipped, tex):
    rows = _supp_names_rows(head, basis)
    if not tex:
        return rows
    return _supp_tex_rows(rows, lambda r: [
        tex_escape(r["partition"]), thousands(r["labels"]),
        thousands(r["accepted"]), thousands(r["synonym"]),
        thousands(r["unresolved"]), r["unresolved_pct"],
        thousands(r["override"]),
        thousands(r["merged"]), thousands(r["named_species"]),
        thousands(r["gbif_added"])], 10)


@supp_table("supp_molecular", ["partition", "category", "f_basis",
                               "source_file", "tips", "molecular_tips",
                               "molecular_fraction", "exact_label_tips",
                               "molecular_species", "product", "f_low",
                               "f_high", "molecular_only_tree",
                               "molecular_only_dated"])
def _t_molecular(head, basis, shipped, tex):
    rows = _supp_molecular_rows(head, basis)
    if not tex:
        return rows

    def render(r):
        src = (r"\texttt{%s}" % tex_escape(r["source_file"])
               if r["source_file"] else "---")
        mt = thousands(r["molecular_tips"])
        if r["molecular_only_tree"] == "yes":
            mt += (r"$^{d}$" if r["molecular_only_dated"] == "yes"
                   else r"$^{g}$")
        return [tex_escape(r["partition"]),
                SUPP_MOLSOURCE_CODE[r["f_basis"]], src, thousands(r["tips"]),
                mt, r["molecular_fraction"],
                _dash(thousands(r["exact_label_tips"])
                      if r["exact_label_tips"] != "" else ""),
                thousands(r["molecular_species"]), r["product"]]
    return _supp_tex_rows(rows, render, 9)


@supp_table("supp_time_scale", ["partition", "category", "roottotip_cv",
                                "ultrametric", "root_height", "time_scale",
                                "study_dating"])
def _t_time(head, basis, shipped, tex):
    rows = _supp_time_rows(head, basis)
    if not tex:
        return rows

    def render(r):
        dates = (TABLE1_DATES_CODE[r["time_scale"]]
                 if r["time_scale"] in ("absolute", "relative")
                 else TABLE1_DATES_CODE[r["study_dating"]])
        return [tex_escape(r["partition"]), r["roottotip_cv"],
                "yes" if r["ultrametric"] == "yes" else "no",
                r["root_height"], dates]
    return _supp_tex_rows(rows, render, 5)


@supp_table("supp_crown", ["partition", "category", "crown_clade",
                           "crown_age_ma", "age_statistic", "interval_lo_ma",
                           "interval_hi_ma", "interval_width_myr",
                           "interval_kind", "n_trees", "calibration",
                           "prior_lo_ma", "prior_hi_ma",
                           "prior_containment_pct", "prior_width_ratio",
                           "prior_q_lo", "prior_q_center", "prior_q_hi",
                           "parent_calibration", "child_calibration",
                           "prior_range", "status"])
def _t_crown(head, basis, shipped, tex):
    rows = _supp_crown_rows(head, basis)
    if not tex:
        return rows

    def render(r):
        iv = ("%s--%s" % (r["interval_lo_ma"], r["interval_hi_ma"])
              if r["interval_kind"] else "---")
        pr = ("%s--%s%s" % (r["prior_lo_ma"], r["prior_hi_ma"],
                            "\\textsuperscript{b}"
                            if r["prior_range"] == "band" else "")
              if r["prior_lo_ma"] else "---")
        return [tex_escape(r["partition"]), tex_escape(r["crown_clade"]),
                r["crown_age_ma"], iv, _dash(r["interval_width_myr"]),
                SUPP_KIND_CODE.get(r["interval_kind"], "---"),
                _dash(thousands(r["n_trees"]) if r["n_trees"] else ""),
                SUPP_CAL_CODE.get(r["calibration"], "---"), pr,
                _dash(r["prior_containment_pct"]),
                _dash(r["prior_width_ratio"]),
                (", ".join(_tex_quantile(r[k]) for k in
                           ("prior_q_lo", "prior_q_center", "prior_q_hi"))
                 if r["prior_q_center"] else "---"),
                SUPP_CAL_CODE.get(r["parent_calibration"], "---"),
                SUPP_CAL_CODE.get(r["child_calibration"], "---"),
                TABLE1_INTERVAL_CODE[r["status"]]]
    return _supp_tex_rows(rows, render, 15)


def _tex_quantile(q):
    """A prior quantile for the Supplement: '<0.025' and '>0.975' as math."""
    if q.startswith("<"):
        return "$<$%s" % q[1:]
    if q.startswith(">"):
        return "$>$%s" % q[1:]
    return "%.2f" % float(q)


@supp_table("supp_sensitivity", ["alternative", "denominator"]
            + list(LEVEL_ORDER))
def _t_sens(head, basis, shipped, tex):
    rows = _supp_sensitivity_rows(head, basis)
    if not tex:
        return rows
    out = []
    for r in rows:
        if r["alternative"] in ("insects_minus", "col_filled", "rep",
                                "drop_automated_synthesis", "bird_band_low",
                                "rule_p_counted"):
            out.append(r"\midrule")
        out.append(" & ".join(
            [SUPP_SENS_ROWS[r["alternative"]], thousands(r["denominator"])]
            + [("$\\leq$ " if r["alternative"] == "union" and k == "C1"
                else "") + _dash(r[k]) for k in LEVEL_ORDER]) + r" \\")
    return "\n".join(out) + "\n"


@supp_table("supp_denominators", ["partition", "category", "described",
                                  "described_source", "col_count",
                                  "col_filled_from", "named_pct",
                                  "named_col_pct", "named_all_trees_pct"])
def _t_den(head, basis, shipped, tex):
    rows = _supp_denominator_rows(head, basis)
    if not tex:
        return rows
    return _supp_tex_rows(rows, lambda r: [
        tex_escape(r["partition"]), thousands(r["described"]),
        tex_escape(r["described_source"]),
        thousands(r["col_count"]) + ("$^{f}$" if r["col_filled_from"] else ""),
        _dash(r["named_pct"]), _dash(r["named_col_pct"]),
        _dash(r["named_all_trees_pct"])], 7)


@supp_table("supp_opentree", ["partition", "category", "ot_phylo_species",
                              "molecular_species", "both", "molecular_not_ot",
                              "both_pct", "molecular_keys", "ot_keys",
                              "keys_both", "keys_both_pct",
                              "keys_both_ot_pct"])
def _t_ot(head, basis, shipped, tex):
    rows = _supp_opentree_rows(head, basis)
    if not tex:
        return rows
    return _supp_tex_rows(rows, lambda r: [
        tex_escape(r["partition"]), thousands(r["ot_phylo_species"]),
        thousands(r["molecular_species"]), thousands(r["both"]),
        thousands(r["molecular_not_ot"]), _dash(r["both_pct"]),
        thousands(r["molecular_keys"]), thousands(r["ot_keys"]),
        thousands(r["keys_both"]), _dash(r["keys_both_pct"]),
        _dash(r["keys_both_ot_pct"])], 11)


SUPP_FIT_LABEL = OrderedDict([
    (("full", "tips"), "full, own $\\rho$"),
    (("molecular", "tips"), "pruned, own $\\rho$"),
    (("molecular", "full_misapplied"), "pruned, full-tree $\\rho$"),
    (("molecular", "family"), "pruned, family $\\rho$"),
])


@supp_table("supp_worked_fits", ["partition", "tree", "rho_convention",
                                 "n_tips", "rho", "lambda", "mu", "r", "r_lo",
                                 "r_hi", "epsilon", "epsilon_lo", "epsilon_hi",
                                 "loglik", "n_events", "n_polytomy_events"])
def _t_fits(head, basis, shipped, tex):
    rows = _supp_fit_rows(head)
    if not tex:
        return rows
    out, prev = [], None
    for r in rows:
        if r["partition"] != prev:
            if prev is not None:
                out.append(r"\midrule")
            prev = r["partition"]
            grp = tex_escape(r["partition"])
        else:
            grp = ""
        out.append(" & ".join([
            grp, SUPP_FIT_LABEL[(r["tree"], r["rho_convention"])],
            thousands(r["n_tips"]), _dash(r["rho"]), r["lambda"], r["mu"],
            "%s (%s--%s)" % (r["r"], r["r_lo"], r["r_hi"]),
            "%s (%s--%s)" % (r["epsilon"], r["epsilon_lo"], r["epsilon_hi"]),
            r["loglik"].replace("-", "$-$"), thousands(r["n_events"]),
            thousands(r["n_polytomy_events"])]) + r" \\")
    return "\n".join(out) + "\n"


SUPP_POLY_OUTCOME = OrderedDict([
    ("rescaled", "rescaled"), ("mu0", "$\\mu = 0$"),
    ("mu0_rejected", "$\\mu = 0$$^{*}$"), ("mixed", "mixed"),
])


@supp_table("supp_worked_polytomy", ["partition", "events_full",
                                     "polytomy_events_full", "events_pruned",
                                     "polytomy_events_pruned", "r_ratio",
                                     "r_ratio_resolved", "epsilon_matching",
                                     "epsilon_matching_resolved", "lr",
                                     "lr_resolved", "outcome_resolved"])
def _t_poly(head, basis, shipped, tex):
    rows = _supp_polytomy_rows(head)
    if not tex:
        return rows
    return _supp_tex_rows(rows, lambda r: [
        tex_escape(r["partition"]), thousands(r["polytomy_events_full"]),
        thousands(r["events_full"]), thousands(r["polytomy_events_pruned"]),
        thousands(r["events_pruned"]), r["r_ratio"], r["r_ratio_resolved"],
        r["epsilon_matching"], r["epsilon_matching_resolved"], r["lr"],
        r["lr_resolved"], SUPP_POLY_OUTCOME[r["outcome_resolved"]]], 12,
        category_key=None)


def _tips_cell(v):
    """Tip counts of a set: one count, or several joined by '/'."""
    if v in ("", None):
        return "---"
    return "/".join(thousands(x) for x in v.split("/"))


def _supp_sets_rows(head):
    """Round F2: one row per archived tree set the fits are repeated on."""
    ts = head["worked_example_S13"]["tree_sets"]

    def rng(a, b, fmt):
        lo, hi = fmt % a, fmt % b
        return lo if lo == hi else "%s--%s" % (lo, hi)
    rows = []
    for v in ts.values():
        rows.append(OrderedDict([
            ("partition", v["partition"]),
            ("set", WORKED_SET_TAG[v["set"]][1]),
            ("n_trees", v["n_trees"]),
            ("tips_pruned", "/".join(str(x) for x in v["n_tips_pruned"])),
            ("tips_full", "/".join(str(x) for x in v["n_tips_full"])),
            ("epsilon_matching", "%s (%s--%s)" % (
                "%.2f" % v["epsilon_matching_median"],
                "%.2f" % v["epsilon_matching_min"],
                "%.2f" % v["epsilon_matching_max"])),
            ("lr", "%s (%s)" % ("%.2f" % v["lr_median"],
                                rng(v["lr_min"], v["lr_max"], "%.2f"))),
            ("n_boundary", v["n_boundary"]),
            ("r_ratio", "" if v["r_ratio_min"] is None else
             rng(v["r_ratio_min"], v["r_ratio_max"], "%.2f")),
            ("summary_epsilon_matching", "%.2f" % v["summary_tree_epsilon_matching"]),
            ("summary_lr", "%.2f" % v["summary_tree_lr"]),
        ]))
    return rows


@supp_table("supp_worked_sets", ["partition", "set", "n_trees", "tips_pruned",
                                 "tips_full", "epsilon_matching", "lr",
                                 "n_boundary", "r_ratio",
                                 "summary_epsilon_matching", "summary_lr"])
def _t_sets(head, basis, shipped, tex):
    rows = _supp_sets_rows(head)
    if not tex:
        return rows
    out, prev = [], None
    for r in rows:
        grp = tex_escape(r["partition"]) if r["partition"] != prev else ""
        prev = r["partition"]
        out.append(" & ".join([
            grp, tex_escape(r["set"]), thousands(r["n_trees"]),
            _tips_cell(r["tips_pruned"]), _tips_cell(r["tips_full"]),
            r["epsilon_matching"], r["lr"], thousands(r["n_boundary"]),
            _dash(r["r_ratio"]), r["summary_epsilon_matching"],
            r["summary_lr"]]) + r" \\")
    return "\n".join(out) + "\n"


@supp_table("supp_worked_families", ["partition", "families",
                                     "tips_with_family", "tips", "f_min",
                                     "f_median", "f_max", "f_pooled",
                                     "rho_uniform", "r_first", "r_last",
                                     "epsilon_first", "epsilon_last"])
def _t_fam(head, basis, shipped, tex):
    rows = _supp_family_rows(head)
    if not tex:
        return rows
    out = []
    for r in rows:
        out.append(" & ".join([
            tex_escape(r["partition"]), thousands(r["families"]),
            thousands(r["tips_with_family"]), thousands(r["tips"]),
            r["f_min"], r["f_median"], r["f_max"], r["f_pooled"],
            r["rho_uniform"], r["epsilon_first"], r["epsilon_last"],
            r["r_first"], r["r_last"]]) + r" \\")
    return "\n".join(out) + "\n"


SUPP_TERMS = [
    # (term in the paper, columns of Data/coverage_basis.csv, other file named
    # by its path in the data release: Supplement S10, release_path_map.csv)
    ("group", ("partition", "category"), ""),
    ("described species", ("described", "described_source"), ""),
    ("reference tree", ("canonical_tree",), "standardized/trees/"),
    ("tips", ("labels_unique",), ""),
    ("named coverage (level C0)", ("c0_units",), ""),
    ("tip labels resolving to no accepted species", ("c0_unresolved_labels",),
     "standardized/col_resolution.csv.gz"),
    ("molecular coverage (level C1)", ("c1_units",), ""),
    ("how the molecular tips are known", ("c1_basis", "f_basis",
                                         "molecular_source_file"), ""),
    ("molecular fraction", ("f_mol", "f_low", "f_high"), ""),
    ("molecular tips", (), "registers/molecular_tips.csv"),
    ("chronogram", ("is_chronogram", "roottotip_cv"), ""),
    ("dated in absolute time (level C1d)", ("time_scale",), ""),
    ("interval coverage (level C1b)", ("interval_recoverable",
                                       "interval_kind"),
     "registers/crown_intervals.csv"),
    ("how the tree was assembled", ("topology_assembly",), ""),
    ("crown age", (), "paper/crown_ages.csv"),
]


@supp_table("supp_terms", ["term", "columns", "file"])
def _t_terms(head, basis, shipped, tex):
    hdr = set(basis[0].keys())
    rows = []
    for t, cols, f in SUPP_TERMS:
        for c in cols:
            if c not in hdr:
                die("Supplement terms table: %s is not a column of %s"
                    % (c, rel(OUT["coverage_basis.csv"])))
        rows.append(OrderedDict([("term", t), ("columns", "; ".join(cols)),
                                 ("file", f)]))
    if not tex:
        return rows
    return "\n".join(r"%s & %s & %s \\" % (
        tex_escape(r["term"]),
        r"\texttt{%s}" % tex_escape(r["columns"]) if r["columns"] else "---",
        r"\texttt{%s}" % tex_escape(r["file"]) if r["file"] else "---")
        for r in rows) + "\n"


@supp_table("supp_checks", ["check", "title"])
def _t_checks(head, basis, shipped, tex):
    rows = [OrderedDict([("check", c), ("title", t)])
            for c, t in _verify_check_list()]
    if not tex:
        return rows
    return "\n".join(r"\texttt{%s} & %s \\" % (tex_escape(r["check"]),
                                               tex_escape(r["title"]))
                     for r in rows) + "\n"


SUPP_OUT_NAMES = ["supp_macros.tex"] + ["%s%s" % (n, ext) for n in SUPP_TABLES
                                        for ext in (".csv", ".tex")]
for _n in SUPP_OUT_NAMES:
    OUT[_n] = os.path.join(SUPP_DIR, _n)
OUT.move_to_end("headline_numbers.json")


def build_headline_json(head):
    return json.dumps(head, indent=2, ensure_ascii=False, sort_keys=False,
                      default=str) + "\n"


# ==========================================================================
# RETIRED FIGURES and OPEN FLAGS
# ==========================================================================

RETIRED = [
    "THE 23.7 % COVERAGE HEADLINE. A hard-coded literal at "
    "Paper/code/figures/fig2_sensitivity_bounds.py:61, reachable from no "
    "numerator/denominator pair in the deposit; its implied denominator, "
    "2,084,409, exists in no file. WITHDRAWN, NOT ADJUSTED. This script owns "
    "no file under Paper/code/, so the deletion of that file and its literal "
    "is REPORTED here, not performed.",

    "THE OLD MIDDLE RUNG 'C1 = 8.65 %'. It never had a principled definition: "
    "the published table applied a molecular fraction to three partitions at "
    "one criterion and to every partially-imputed canonical at the next, a "
    "scope not derivable from any deposited field.",

    "THE OLD STRICT CRITERION, which dropped whole partitions categorically "
    "while scaling others proportionally -- two different operations under one "
    "label. It could not tell 85 % molecular support from 5 %: both scored "
    "zero. The reader who wants that view now reads it off ONE sensitivity "
    "line (sensitivities.leave_assembly_class_out.retain_only_primary_analysis).",

    "THE min(x, 100) COVERAGE CAP and the hard-coded collection-backed slug "
    "list in prep_data.py. No eukaryotic canonical partition exceeds 100 %, "
    "proved by the invariant check on every run, so a cap could only ever mask "
    "a defect. This script owns neither file; the deletions are reported.",

    "THE LEDGER'S `dated` COLUMN as a source of truth. Time calibration is now "
    "measured from tree geometry and shipped with its root-to-tip CV, so a "
    "reader may set their own tolerance instead of inheriting a bit.",
]

OPEN_FLAGS = [
    "PER-TIP ENUMERATION IS PARTIAL. molecular_tips.csv enumerates C1 only for "
    "the partitions whose molecular source is a deposited tree. Everywhere else "
    "C1 is an expectation, so the paper must write TIP-EQUIVALENTS, never "
    "species, at C1 and below. Completing the enumeration is the single "
    "highest-value remaining data task.",

    "LABEL-MATCHING RULE FOR f_mol -- A DELIBERATE CHOICE, FLAGGED FOR REVIEW. "
    "The intersection behind f_mol resolves BOTH sides through the deposit's "
    "own name mapping (dictionary.csv for the shipped tree's numeric ids; "
    "full_mapping.csv for the raw source's original labels). Naive exact-label "
    "matching was rejected because it is a measurement artefact: it scores a "
    "source whose labels carry collector or voucher decorations at or near "
    "ZERO, which is a fact about label formatting, not about molecular support. "
    "The adopted rule changes exactly one partition's number relative to naive "
    "matching, and it moves that number CLOSER to the numerator the ground "
    "truth independently recorded. Both figures are published per partition in "
    "headline_numbers.json -> molecular_evidence "
    "(measured_intersection vs exact_label_intersection), so either can be "
    "reproduced. If the exact-label rule is preferred instead, that is a "
    "one-line change here and it must be made consistently for every "
    "partition, not only where it is convenient.",

    "TOPOLOGY_ASSEMBLY RESTS ON UNRECONCILED LEDGER FIELDS. It is derived from "
    "`inference` and `otl_pipeline_involvement` in Data/tree_provenance.csv -- "
    "the same ledger whose `dated` column is wrong wherever it can be checked. "
    "That is exactly why it is a NOMINAL ANNOTATION and never a rung, and "
    "Methods must say so. If those two fields are ever reconciled and "
    "certified, the annotation can be promoted; until then it cannot.",

    "TWO TREE COLLECTIONS. Data/site/data/partitions/<slug>/tree.nwk is a "
    "second, curated collection that disagrees with Data/standardized/trees/ on "
    "several partitions. The truth files and therefore this script follow the "
    "standardized collection, per the settled 'the shipped file governs' rule. "
    "If the curated collection is the deposit of record instead, every number "
    "here must be rebuilt.",

    "FIGURE-1 WORKAROUND NOW REDUNDANT. The deposited fig1_coverage_bars.csv "
    "carried several canonical eukaryote partitions with is_unrepresented=True "
    "and blank tip counts, and Paper/code/figures/fig1_dotplot.py compensates "
    "by overriding them from a separate overlay CSV. The rebuilt file fixes the "
    "root cause, so that reconciliation block should be DELETED rather than "
    "left to fight a file that is now correct.",

    "SUPERSEDED WORKING COPIES LEFT IN PLACE. "
    "Paper/figure2-3-rework-2026-05-31/c0c1c2_per_partition.csv belongs to a "
    "dated working session and to another owner; it was not touched, but it now "
    "holds superseded numbers under a retired criterion. Any figure code still "
    "reading it must be repointed at Paper/inputs/coverage_per_partition.csv.",

    "NEMATODES ATTRIBUTION. The main table prints the corrected count but "
    "credits the Catalogue of Life, whose figure is the rejected one. The value "
    "used here carries truth_species.csv's own attribution. The fix is owed in "
    "a .tex file this script does not own.",
]


# ==========================================================================
# DIFF REPORTING
# ==========================================================================

def summarise_change(path, new_text, quiet):
    """A compact, human-readable before/after summary for one file."""
    name = rel(path)
    if not os.path.exists(path):
        print("  NEW      %s  (%d lines)" % (name, new_text.count("\n")))
        return
    with open(path, encoding="utf-8", errors="replace") as fh:
        old = fh.read()
    if old == new_text:
        print("  same     %s" % name)
        return
    o, n = old.count("\n"), new_text.count("\n")
    print("  CHANGED  %s  (%d -> %d lines)" % (name, o, n))
    if quiet:
        return
    if name.endswith(".json"):
        for line in _diff_json(old, new_text):
            print("             " + line)
    elif name.endswith(".csv"):
        for line in _diff_csv_cells(old, new_text):
            print("             " + line)
    else:
        shown = 0
        for line in difflib.unified_diff(old.splitlines(), new_text.splitlines(),
                                         lineterm="", n=0):
            if line.startswith(("---", "+++", "@@")):
                continue
            print("             " + line[:160])
            shown += 1
            if shown >= 12:
                print("             ... (further differences suppressed)")
                break


def _diff_csv_cells(old_text, new_text, max_show=25):
    """Cell-level diff keyed on the first column, so a reader sees WHICH
    partition changed and by how much, not a wall of reordered lines."""
    def load(t):
        """Key each row by the fewest leading columns that make it unique.

        Keying on the first column alone silently collapses tables whose first
        column repeats (a category x type tally, for instance), and then reports
        one row's change as though it were another's.
        """
        rows = list(csv.reader(io.StringIO(t)))
        if not rows:
            return [], OrderedDict()
        hdr, body = rows[0], [r for r in rows[1:] if r]
        for width in range(1, max(len(hdr), 1) + 1):
            keys = ["|".join(r[:width]) for r in body]
            if len(set(keys)) == len(keys):
                return hdr, OrderedDict(zip(keys, body))
        # Not uniquely keyable on any prefix: fall back to position.
        return hdr, OrderedDict(("row %d|%s" % (i, "|".join(r[:2])), r)
                                for i, r in enumerate(body, 1))

    oh, orows = load(old_text)
    nh, nrows = load(new_text)
    out = []
    if oh != nh:
        added = [c for c in nh if c not in oh]
        removed = [c for c in oh if c not in nh]
        if added:
            out.append("columns added:   %s" % ", ".join(added))
        if removed:
            out.append("columns removed: %s" % ", ".join(removed))
    gone = [k for k in orows if k not in nrows]
    new = [k for k in nrows if k not in orows]
    if gone:
        out.append("rows removed (%d): %s" % (len(gone), ", ".join(gone[:8])))
    if new:
        out.append("rows added (%d): %s" % (len(new), ", ".join(new[:8])))
    shared_cols = [c for c in nh if c in oh]
    changes = 0
    for key in nrows:
        if key not in orows:
            continue
        orow = dict(zip(oh, orows[key]))
        nrow = dict(zip(nh, nrows[key]))
        diffs = ["%s %s -> %s" % (c, orow.get(c, ""), nrow.get(c, ""))
                 for c in shared_cols if orow.get(c, "") != nrow.get(c, "")]
        if diffs:
            changes += 1
            if changes <= max_show:
                out.append("%s: %s" % (key, "; ".join(diffs)[:190]))
    if changes > max_show:
        out.append("... %d further changed row(s) suppressed"
                   % (changes - max_show))
    if not out:
        out.append("(row order or whitespace only)")
    return out


def _flatten_json(obj, prefix=""):
    flat = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            flat.update(_flatten_json(v, "%s.%s" % (prefix, k) if prefix else k))
    elif isinstance(obj, list):
        flat[prefix] = "[%d item(s)]" % len(obj)
    else:
        flat[prefix] = obj
    return flat


def _diff_json(old_text, new_text, max_show=30):
    try:
        o = _flatten_json(json.loads(old_text))
        n = _flatten_json(json.loads(new_text))
    except ValueError:
        return ["(previous file was not valid JSON; replaced wholesale)"]
    out, shown = [], 0
    for k in sorted(set(o) | set(n)):
        if k.startswith("warnings") or k.endswith("notes"):
            continue
        if o.get(k) != n.get(k):
            shown += 1
            if shown <= max_show:
                out.append("%s: %r -> %r" % (k, o.get(k), n.get(k)))
    if shown > max_show:
        out.append("... %d further key(s) changed" % (shown - max_show))
    return out or ["(no scalar key changed)"]


# ==========================================================================
# INVARIANTS -- structural checks, not value assertions.  No expected number
# is written down here; each check is a property that must hold whatever the
# data say.
# ==========================================================================

def check_invariants(truth, basis, head, shipped, warnings):
    problems = []
    euk = [r for r in basis if is_euk(r)]
    den = head["denominators"]["eukaryotic"]["described"]

    # 1. The basis table is the whole partition set, once each.
    if len(basis) != len(truth.species):
        problems.append(
            "coverage_basis.csv has %d rows for %d partitions; the basis table "
            "must carry EVERY partition, including the treeless ones, or a "
            "denominator can be formed that is not a column sum of it."
            % (len(basis), len(truth.species)))
    if len({r["partition"] for r in basis}) != len(basis):
        problems.append("coverage_basis.csv repeats a partition.")

    # 2. Monotonicity, PER PARTITION, over the nested weights.
    order = list(RUNGS)
    for r in euk:
        ws = [rung_weight(r, k) for k in order]
        for a, b in zip(ws, ws[1:]):
            if b > a + 1e-12:
                problems.append(
                    "partition %r violates rung monotonicity (weights %s): a "
                    "stricter rung cannot admit more than a looser one."
                    % (r["partition"], [round(w, 6) for w in ws]))
                break

    # 3. Monotonicity of the totals.
    tot = [head["rungs"][k]["numerator_sum"] for k in order]
    for (ka, a), (kb, b) in zip(list(zip(order, tot)), list(zip(order, tot))[1:]):
        if b > a + 1e-6:
            problems.append("rung totals are not monotone: %s %.3f < %s %.3f"
                            % (ka, a, kb, b))

    # 4. No partition exceeds 100 % -- which is WHY the cap was deleted.
    for r in euk:
        d = as_int(r["described"], 0)
        t = as_int(r["c0_units"], 0)
        if d and t > d:
            problems.append(
                "partition %r has %d distinct tip labels against %d described "
                "species (>100 %% coverage). The min(x,100) cap was deleted "
                "because no partition did this; one now does, and capping it "
                "would MASK the defect."
                % (r["partition"], t, d))

    # 5. C0's union can never exceed its sum, nor the denominator.
    c0 = head["rungs"]["C0"]
    if c0["numerator_union"] > c0["numerator_sum_rounded"]:
        problems.append("the C0 union (%d) exceeds the per-partition sum (%d)."
                        % (c0["numerator_union"], c0["numerator_sum_rounded"]))
    if c0["numerator_sum_rounded"] > den:
        problems.append("C0 (%d) exceeds the eukaryotic denominator (%d)."
                        % (c0["numerator_sum_rounded"], den))

    # 6. Denominator identity: it must BE the column sum of the basis file.
    if sum(as_int(r["described"], 0) for r in euk) != den:
        problems.append(
            "the reported eukaryotic denominator is not the `described` column "
            "sum of coverage_basis.csv. This is THE structural guarantee "
            "against a phantom denominator and it has failed.")

    # 7. No prokaryote may appear in a eukaryotic numerator.
    for r in basis:
        if not is_euk(r) and rung_weight(r, "C0") > 0 and r in euk:
            problems.append("a prokaryotic partition entered the eukaryotic "
                            "numerator: %r" % r["partition"])

    # 8. Every rung percentage must be reproducible from its own numerator and
    #    denominator -- i.e. no percentage anywhere is an independent literal.
    for k in order:
        e = head["rungs"][k]
        want = pct(e["numerator_sum"], e["denominator"])
        if want is None or abs(want - e["percent"]) > 1e-9:
            problems.append("rung %s's percentage is not numerator/denominator."
                            % k)

    # 9. f_mol must lie in [0,1] and inside its own band.
    for r in euk:
        if as_int(r["c0_units"], 0) <= 0:
            continue
        f = as_float(r["f_mol"])
        lo, hi = as_float(r["f_low"]), as_float(r["f_high"])
        if f is None or not (0.0 <= f <= 1.0):
            problems.append("partition %r has f_mol=%r outside [0,1]."
                            % (r["partition"], r["f_mol"]))
        elif lo is None or hi is None or not (lo <= f <= hi):
            problems.append("partition %r has f_mol=%s outside its band [%s, %s]."
                            % (r["partition"], r["f_mol"], r["f_low"], r["f_high"]))

    # 10. A treeless partition must contribute nothing to any rung.
    for r in euk:
        if as_int(r["c0_units"], 0) == 0:
            if any(rung_weight(r, k) > 0 for k in order):
                problems.append("treeless partition %r carries a non-zero rung "
                                "weight." % r["partition"])

    # 11. C1b must be a subset condition of C1d: an interval is only meaningful
    #     on a chronogram.
    for r in euk:
        if r.get("interval_recoverable") == "yes" and r.get("is_chronogram") != "yes":
            problems.append(
                "partition %r is marked interval_recoverable but its canonical "
                "tree is not a chronogram: the interval cannot be attached to "
                "the tree the atlas distributes, so archival_class is at the "
                "wrong grain." % r["partition"])

    return problems


# ==========================================================================
# MAIN
# ==========================================================================

def main(argv=None):
    global PARSE_NEWICK

    ap = argparse.ArgumentParser(
        description="Regenerate every derived atlas artefact from "
                    "Data/_rebuild/truth_*.csv and the shipped trees.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="compute everything and print the diff, write nothing")
    ap.add_argument("--only", nargs="+", metavar="NAME", default=None,
                    help="write only these outputs (short names)")
    ap.add_argument("--quiet", action="store_true",
                    help="one summary line per file; suppress cell-level diffs")
    ap.add_argument("--crown-ages-only", action="store_true",
                    help="write only Paper/inputs/crown_ages.csv, before the "
                         "worked example is re-run (its R script reads the "
                         "crown ages; the full rebuild refuses fits computed "
                         "on other crown ages)")
    ap.add_argument("--no-verify-trees", action="store_true",
                    help="skip re-measuring root-to-tip CV against truth_dated.csv")
    args = ap.parse_args(argv)

    if args.only:
        unknown = [n for n in args.only if n not in OUT]
        if unknown:
            die("unknown --only target(s): %s\n                known: %s"
                % (", ".join(unknown), ", ".join(OUT)))

    print("=" * 78)
    print("phylo-atlas rebuild -- the coverage ladder over one basis table"
          + ("  [DRY RUN]" if args.dry_run else ""))
    print("repository: %s" % ROOT)
    print("=" * 78)

    warnings = []
    PARSE_NEWICK = load_audit_parser()
    print("\nparser: imported parse_newick from %s" % rel(AUDIT_COMMON))

    # ---- measure, then load ----------------------------------------------
    print("\nmeasured truth (re-measured from the shipped files)")
    measured = refresh_measured_truth(args.dry_run, warnings)
    truth = Truth(measured)
    print("\ntruth files")
    attr = {"dated": "dated_rows", "tips": "tips_rows",
            "species": "species_rows", "provenance": "prov_rows"}
    for name, path in TRUTH_FILES.items():
        print("  %-11s %-42s %4d rows"
              % (name, rel(path), len(getattr(truth, attr[name]))))
    print("  %d shipped trees | %d partitions (%d eukaryotic, %d prokaryotic) | "
          "%d canonical (%d eukaryotic) | %d alternative"
          % (len(truth.dated), len(truth.species), len(truth.eukaryotic),
             len(truth.species) - len(truth.eukaryotic), len(truth.canon_prov),
             len(truth.canon_euk), len(truth.alt_prov)))

    names = Names(warnings)
    shipped = ShippedTrees(truth, names, warnings, not args.no_verify_trees)
    if not args.no_verify_trees:
        checked = sum(1 for v in shipped.cv_measured.values() if v is not None)
        print("  re-measured root-to-tip CV on %d canonical trees "
              "(tolerance %.0e)" % (checked, CV_RECHECK_TOL))

    ledger = PreAuditLedger(truth, warnings)
    pt = Passthrough(truth, warnings)
    if pt._yaml is None:
        die("PyYAML is not importable. Without it estimated_total and "
            "confidence silently fall back to data_estimates.csv and the "
            "outputs differ from the committed ones. Run with an interpreter "
            "that has PyYAML (pip install pyyaml).")

    print("\nmolecular evidence (axis A, measured)")
    mol = resolve_molecular(truth, shipped, names, warnings)
    for part in sorted(mol):
        if truth.is_prokaryotic(part):
            continue
        ev = mol[part]
        band = ("  band [%.4f, %.4f]" % (ev.f_low, ev.f_high)
                if ev.f_high > ev.f_low else "")
        src = ("  <- %s" % ev.source_file) if ev.source_file else ""
        print("  %-14s f=%.6f  %-24s%s%s"
              % (part, ev.f_mol, ev.f_basis, src, band))

    # ---- compute ----------------------------------------------------------
    basis = build_basis(truth, shipped, mol, ledger, warnings)
    head = compute_headline(truth, basis, shipped, mol, ledger, pt, warnings)

    problems = check_invariants(truth, basis, head, shipped, warnings)
    if problems:
        sys.stderr.write("\nINVARIANT FAILURES -- NOTHING WAS WRITTEN:\n")
        for pb in problems:
            sys.stderr.write("  * %s\n" % pb)
        sys.exit(1)

    # ---- report -----------------------------------------------------------
    d = head["denominators"]
    print("\ndenominators (column sums of the basis table, not literals)")
    print("  eukaryotic  %2d partitions  %s described   <- THE denominator"
          % (d["eukaryotic"]["n_partitions"], commas(d["eukaryotic"]["described"])))
    print("  prokaryotic %2d partitions  %s described   <- own line, never mixed"
          % (d["prokaryotic"]["n_partitions"], commas(d["prokaryotic"]["described"])))
    print("  all life    %2d partitions  %s described   <- SI and Figure 1 only"
          % (d["all_life"]["n_partitions"], commas(d["all_life"]["described"])))

    print("\nthe ladder")
    for k, e in head["rungs"].items():
        extra = ""
        if "numerator_union" in e:
            extra = ("  [union %s = %.4f %%; %d labels shared between two "
                     "partitions]" % (commas(e["numerator_union"]),
                                      e["percent_union"], e["shared_labels"]))
        print("  %-4s %14s %-15s %7.4f %%%s"
              % (k, commas(e["numerator_sum_rounded"]), e["unit"],
                 e["percent"], extra))

    c0 = head["rungs"]["C0"]
    dark = head["disclosures"]["dark_matter"]
    print("\ndark matter (formed on the C0 UNION throughout)")
    print("  %s species (%.1f %%) have no tip anywhere in the atlas"
          % (commas(dark["no_tip_anywhere"]), dark["no_tip_anywhere_pct"]))
    print("  %s species (%.1f %%) have no molecularly supported placement"
          % (commas(dark["no_molecular_placement"]),
             dark["no_molecular_placement_pct"]))

    t = head["trees"]
    print("\nstructure")
    print("  %d of %d canonical trees are measured chronograms; everything at "
          "C1d and below rests on those"
          % (t["n_partition_canonical_chronograms"], t["n_partition_canonical"]))
    print("  archival classes: %s"
          % ", ".join("%s %d" % (k, v)
                      for k, v in head["archival_classes_canonical"].items()))
    print("  topology assembly (NOMINAL annotation): %s"
          % ", ".join("%s %d" % (k, v)
                      for k, v in head["topology_assembly"]["counts"].items()))

    print("\nsensitivities")
    for part, s in head["sensitivities"]["f_mol_bands"].items():
        print("  %s f in [%.3f, %.3f]: C1 %.4f %% -> %.4f %%; %s"
              % (part, s["f_low"], s["f_high"],
                 s["at_f_low"]["C1"]["percent"], s["at_f_high"]["C1"]["percent"],
                 "PROPAGATES past C1" if s["propagates_to_C1d"]
                 else "does NOT propagate to C1d/C1b"))
    la = head["sensitivities"]["leave_assembly_class_out"]
    key = "retain_only_" + ASSEMBLY_PRIMARY
    if key in la:
        r = la[key]["rungs"]
        print("  retain only %s: C0 %.2f %%, C1 %.2f %%, C1d %.2f %%, C1b %.3f %%"
              % (ASSEMBLY_PRIMARY, r["C0"]["percent"], r["C1"]["percent"],
                 r["C1d"]["percent"], r["C1b"]["percent"]))
    ds = head["sensitivities"]["denominator"]
    ro = ds["represented_partitions_only"]
    print("  represented partitions only (D = %s): C0 %.2f %%, C1 %.2f %%, "
          "C1d %.2f %%" % (commas(ro["described"]), ro["rungs"]["C0"]["percent"],
                           ro["rungs"]["C1"]["percent"],
                           ro["rungs"]["C1d"]["percent"]))
    lsc = ds["largest_single_count"]
    print("  elasticity: %s alone is %.1f %% of the denominator; +/-%d %% moves "
          "C0 across %s"
          % (lsc["partition"], lsc["share_of_denominator_pct"],
             int(round(100 * ELASTICITY_FRACTION)),
             " / ".join("%.2f %%" % v["rungs"]["C0"]["percent"]
                        for v in ds["elasticity_demonstration"]["results"].values())))

    disc = head["disclosures"]
    print("\ndisclosures the Results MUST print")
    dom = disc["dominance"]
    print("  %s is %.1f %% of C0, %.1f %% of C1, %.1f %% of C1d -- state it and "
          "frame the asymmetry as the finding"
          % (dom["partition"], dom["share_of_C0_pct"], dom["share_of_C1_pct"],
             dom["share_of_C1d_pct"]))
    for b, v in disc["f_basis_accounting_of_C1"].items():
        print("    C1 by basis: %-24s %9s (%.1f %%)"
              % (b, commas(v["tip_equivalents"]), v["share_of_C1_pct"]))
    sz = disc["structural_zeros"]
    print("  %d of %d eukaryotic partitions (%s species, %.1f %%) have NO TREE "
          "at all -- a structural zero"
          % (sz["n_partitions"], sz["of_eukaryotic_partitions"],
             commas(sz["described"]), sz["share_of_denominator_pct"]))
    ov = disc["overlapping_labels"]
    print("  %d labels occur in two eukaryotic canonical partitions each "
          "(%.3f %% of C0) -- outgroups and homonyms"
          % (ov["n_labels_in_two_or_more_partitions"],
             ov["share_of_C0_union_pct"]))

    # ---- render -----------------------------------------------------------
    # headline_numbers.json renders LAST so it carries the complete warning
    # list; every other builder may append to `warnings`.
    builders = OrderedDict([
        ("coverage_basis.csv", lambda: build_coverage_basis(basis, warnings)),
        ("molecular_tips.csv",
         lambda: build_molecular_tips(truth, mol, shipped, warnings)),
        ("atlas_coverage.py", lambda: ATLAS_COVERAGE_PY),
        ("atlas_coverage.R", lambda: ATLAS_COVERAGE_R),
        ("partition_summary.csv",
         lambda: build_partition_summary(truth, pt, warnings)),
        ("coverage_per_partition.csv",
         lambda: build_coverage_per_partition(basis, truth, warnings)),
        ("fig1_coverage_bars.csv",
         lambda: build_fig1_coverage_bars(truth, pt, warnings)),
        ("fig2_coverage_scatter.csv", lambda: build_fig2_scatter(basis, warnings)),
        ("fig3_c0c1_bars.csv", lambda: build_fig3_c0c1_bars(basis, warnings)),
        ("fig3_dated_undated_counts.csv", lambda: build_fig3_counts(truth, warnings)),
        ("table_s5_tree_provenance.csv",
         lambda: build_table_s5_csv(truth, basis, mol, ledger, warnings)),
        ("table_s5_tree_provenance.tex",
         lambda: build_table_s5_tex(truth, basis, mol, ledger, warnings)),
        ("table_s6_figure1_source_data.csv",
         lambda: build_table_s6_fig1_csv(truth, pt, warnings)),
        ("table_s6_figure1_source_data.tex",
         lambda: build_table_s6_fig1_tex(truth, pt, warnings)),
        ("table_s6_alternative_trees.csv",
         lambda: build_table_s6_alternative_trees(truth, warnings)),
        ("table_s7_richness_estimates.csv",
         lambda: build_table_s7_csv(head, warnings)),
        ("table_s7_richness_estimates.tex",
         lambda: build_table_s7_tex(head, warnings)),
        ("c0_name_resolution.csv",
         lambda: build_c0_name_resolution(head, warnings)),
        ("crown_ages.csv", lambda: build_crown_ages_csv(head, warnings)),
        ("crown_ages_macros.tex",
         lambda: build_crown_ages_macros(head, warnings)),
        ("opentree_comparison.csv", lambda: build_opentree_csv(head, warnings)),
        ("table_datasets.csv",
         lambda: build_dataset_table_csv(truth, basis, warnings)),
        ("table_datasets.tex",
         lambda: build_dataset_table_tex(truth, basis, warnings)),
        ("table_datasets_macros.tex",
         lambda: build_dataset_table_macros(truth, basis, warnings)),
        ("numbers_macros.tex",
         lambda: build_number_macros(head, basis, warnings)),
        ("table_which_tree.csv",
         lambda: build_which_tree_csv(head, basis, warnings)),
        ("table_which_tree.tex",
         lambda: build_which_tree_tex(head, basis, warnings)),
        ("table_rungs_sensitivity.csv",
         lambda: build_rungs_sensitivity_csv(head, warnings)),
        ("table_rungs_sensitivity.tex",
         lambda: build_rungs_sensitivity_tex(head, warnings)),
        ("table_worked_example.csv",
         lambda: build_worked_table_csv(head, warnings)),
        ("table_worked_example.tex",
         lambda: build_worked_table_tex(head, warnings)),
        ("fig1_coverage_cascade.csv",
         lambda: build_fig_cascade_csv(head, basis, warnings)),
        ("supp_macros.tex",
         lambda: build_supp_macros(head, basis, shipped, warnings)),
    ] + [("%s%s" % (n, ext),
          (lambda n=n, ext=ext: (build_supp_csv if ext == ".csv"
                                 else build_supp_tex)(n, head, basis, shipped)))
         for n in SUPP_TABLES for ext in (".csv", ".tex")] + [
        ("headline_numbers.json", None),      # rendered below
    ])

    # EVERY builder always runs, even when --only narrows what is written:
    # builders append to `warnings`, and headline_numbers.json embeds that list,
    # so rendering a subset would make the JSON depend on which targets were
    # requested and break idempotency.  Only the WRITE step is narrowed.
    head["revision_checks"] = revision_checks(truth, shipped, names, warnings)
    head["name_resolution_S05"] = name_resolution(truth, shipped, basis)
    head["crown_ages"] = crown_ages(truth, basis, warnings)
    head["coverage_union"] = union_level_assignment(head, basis, shipped)
    head["union_bound_S09"] = union_bound(truth, basis, shipped, warnings)
    head["opentree_P01"] = opentree_comparison(basis, shipped, mol)
    if args.crown_ages_only:
        # the inputs of worked_example_bd.R that the rebuild writes: the crown
        # ages and the basis table (described species, tips)
        for name, txt in (("crown_ages.csv", build_crown_ages_csv(head, warnings)),
                          ("coverage_basis.csv",
                           build_coverage_basis(basis, warnings))):
            summarise_change(OUT[name], txt, args.quiet)
            if not args.dry_run:
                with open(OUT[name], "w", encoding="utf-8") as fh:
                    fh.write(txt)
        if not args.dry_run:
            print("wrote %s and %s only; now run Rscript "
                  "Paper/code/analysis/worked_example_bd.R and the full "
                  "rebuild" % (rel(OUT["crown_ages.csv"]),
                               rel(OUT["coverage_basis.csv"])))
        return 0
    head["worked_example_S13"] = worked_example(basis, mol, shipped,
                                                head["crown_ages"], warnings)
    rendered = OrderedDict()
    for name, fn in builders.items():
        if fn is not None:
            rendered[name] = fn()
    head["warnings"] = warnings
    rendered["headline_numbers.json"] = build_headline_json(head)

    # ---- diff + write -----------------------------------------------------
    targets = args.only or list(OUT)
    print("\nbefore / after")
    for name in targets:
        summarise_change(OUT[name], rendered[name], args.quiet)

    print("\nwriting" + (" (dry run)" if args.dry_run else ""))
    log, results = [], defaultdict(list)
    for name in targets:
        results[atomic_write(OUT[name], rendered[name], args.dry_run,
                             log, warnings)].append(rel(OUT[name]))
    for line in log:
        print(line)
    for status in sorted(results):
        print("  %-13s %d file(s)" % (status + ":", len(results[status])))

    print("\nnot rebuilt (left exactly as deposited)")
    for path, why in NOT_REBUILT.items():
        print("  %-44s %s" % (path, why))

    if warnings:
        print("\nwarnings (%d)" % len(warnings))
        for w in warnings:
            print("  ! %s" % w)

    print("\nRETIRED, not adjusted (%d)" % len(RETIRED))
    for i, r in enumerate(RETIRED, 1):
        print("  %2d. %s" % (i, re.split(r"\.\s", r, maxsplit=1)[0][:180]))

    print("\nSTILL OPEN -- flagged, not decided (%d)" % len(OPEN_FLAGS))
    for i, f in enumerate(OPEN_FLAGS, 1):
        print("  %2d. %s" % (i, re.split(r"\.\s", f, maxsplit=1)[0][:180]))
    print("\n  full text of every flag is in %s under `open_flags`."
          % rel(OUT["headline_numbers.json"]))

    print("\ndone." + (" (dry run: no file was modified)" if args.dry_run else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
