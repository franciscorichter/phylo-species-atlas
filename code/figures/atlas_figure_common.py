"""Shared loaders, path resolution and palette for the atlas figure scripts.

Design rules this module exists to enforce (all four are audit findings that
this rebuild is meant to make structurally impossible to reintroduce):

1.  NO ABSOLUTE PATHS.  Every root is derived from ``__file__``.  The scripts
    therefore run from any checkout, on any machine, from any working
    directory, and the two synced copies of the tree behave identically
    without a local symlink.

2.  NO HARD-CODED QUANTITIES.  Every number a figure draws or prints is read
    from the deposit:
      * ``Data/coverage_basis.csv``          -- the per-partition basis; the
        column sums of this file ARE the numerators and denominators.
      * ``Paper/inputs/headline_numbers.json`` -- the generated headline file,
        used here only to CROSS-CHECK what we recompute from the basis.
    A literal such as the withdrawn 23.7 % headline cannot survive, because
    there is nowhere left to write one down.

3.  THE LADDER IS THE CURRENT ONE.  C0 NAMED / C1 MOLECULAR / C1d DATED /
    C1b BOUNDED.  The old permissive/molecular/strict "C2" triple is retired:
    it dropped whole partitions categorically while scaling others
    proportionally, so it could not tell 85 % molecular support from 5 %.

4.  TIME CALIBRATION IS MEASURED, NOT INHERITED.  ``is_chronogram`` in the
    basis file is the root-to-tip-CV measurement on the shipped ``.nwk`` in
    absolute time (decision D14), not the ledger's ``dated`` bit; the count
    is ``trees.n_partition_canonical_chronograms`` in the headline file.

Every consumer should call :func:`check_against_headlines` after computing its
own aggregates.  It raises rather than warns: a figure that silently disagrees
with the deposit is worse than a figure that fails to build.
"""
from __future__ import annotations

import csv
import json
import math
import pathlib
import re

# ── Roots, derived from this file's location ────────────────────────────────
# .../phylo-atlas/Paper/code/figures/atlas_figure_common.py
#                 ^^^^^ PAPER_ROOT is three parents up.
FIGURES_CODE_DIR = pathlib.Path(__file__).resolve().parent
PAPER_ROOT = FIGURES_CODE_DIR.parent.parent
REPO_ROOT = PAPER_ROOT.parent
DATA_ROOT = REPO_ROOT / "Data"
INPUTS_DIR = PAPER_ROOT / "inputs"
FIGURES_DIR = PAPER_ROOT / "figures"

# In the public data repository (github.com/franciscorichter/phylo-species-atlas)
# this file sits at code/figures/: the generated inputs are in paper/, the
# per-group table in registers/, and figures are written to paper/figures/.
if (FIGURES_CODE_DIR.parent.parent / "registers").is_dir():
    REPO_ROOT = FIGURES_CODE_DIR.parent.parent
    PAPER_ROOT = REPO_ROOT / "paper"
    DATA_ROOT = REPO_ROOT / "registers"
    INPUTS_DIR = PAPER_ROOT
    FIGURES_DIR = PAPER_ROOT / "figures"

COVERAGE_BASIS = DATA_ROOT / "coverage_basis.csv"
HEADLINE_JSON = INPUTS_DIR / "headline_numbers.json"

RUNGS = ("C0", "C1", "C1d", "C1b")

# ── Palette ────────────────────────────────────────────────────────────────
# Okabe-Ito, colour-blind safe.  One palette for the whole figure set so a
# rung keeps its colour from one panel to the next.
RUNG_COLOR = {
    "C0":  "#0072B2",   # blue     -- NAMED
    "C1":  "#E69F00",   # orange   -- MOLECULAR
    "C1d": "#009E73",   # green    -- DATED
    "C1b": "#CC79A7",   # magenta  -- BOUNDED
}
DESCRIBED_COLOR = "#333333"
INTERVAL_COLOR = "#7f9dbf"
BAND_COLOR = "#999999"
PROKARYOTE_COLOR = "#7d7d7d"

RUNG_SHORT = {
    "C0":  "NAMED",
    "C1":  "MOLECULAR",
    "C1d": "DATED",
    "C1b": "BOUNDED",
}


def apply_style(mpl, size=8):
    """House rcParams.  Helvetica with a real fallback chain."""
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = [
        "Helvetica", "Helvetica Neue", "Arial", "DejaVu Sans",
    ]
    mpl.rcParams["font.size"] = size
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["pdf.fonttype"] = 42          # editable text in the PDF
    mpl.rcParams["ps.fonttype"] = 42


# ── Small parsers ──────────────────────────────────────────────────────────
def _int(x, default=0):
    x = (x or "").strip()
    if not x:
        return default
    try:
        return int(float(x))
    except ValueError:
        return default


def _float(x, default=None):
    x = (x or "").strip()
    if not x:
        return default
    try:
        return float(x)
    except ValueError:
        return default


def _yes(x):
    return (x or "").strip().lower() in ("yes", "true", "1")


def slugify(name: str) -> str:
    """Partition display name -> the slug convention used across the deposit."""
    return (name or "").strip().lower().replace(" ", "_").replace("&", "and")


# ── The basis file ─────────────────────────────────────────────────────────
def load_basis(path: pathlib.Path | None = None) -> list[dict]:
    """Load Data/coverage_basis.csv, typed, with the four rungs derived.

    One row per partition (62: 60 eukaryotic + 2 prokaryotic).  The derived
    fields are computed here, once, so that no figure invents its own rule:

        C0  = c0_units                           (accepted names; NAMED)
        C1  = c1_units: distinct accepted names among molecular tips
        C1d = C1 if the shipped tree is a MEASURED chronogram else 0
        C1b = C1d if the divergence-time interval is recoverable else 0

    ``f_low``/``f_high`` carry the molecular-fraction band where the deposit
    records one; where it does not they equal ``f_mol``.
    """
    path = path or COVERAGE_BASIS
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            described = _int(r.get("described"))
            labels = _int(r.get("labels_unique"))
            # C0 counts distinct Catalogue of Life accepted names (S05, D18)
            # for eukaryotes and GTDB species clusters for prokaryotes (D12);
            # the rebuild writes that count to c0_units.
            c0_units = _int(r.get("c0_units"))
            f_mol = _float(r.get("f_mol"))
            has_tree = bool((r.get("canonical_tree") or "").strip())
            chrono = _yes(r.get("is_chronogram"))
            recoverable = _yes(r.get("interval_recoverable"))

            c0 = float(c0_units) if has_tree else 0.0
            # C1 counts distinct accepted names among the molecular tips
            # (c1_units, audit A01); c0 * f_mol only where the basis has none
            c1_units = _float(r.get("c1_units"))
            c1 = (c1_units if (has_tree and c1_units is not None)
                  else c0 * (f_mol if f_mol is not None else 0.0))
            c1d = c1 if chrono else 0.0
            c1b = c1d if recoverable else 0.0

            rows.append({
                "partition": (r.get("partition") or "").strip(),
                "group": (r.get("group") or "").strip(),
                "slug": slugify(r.get("partition")),
                "category": (r.get("category") or "").strip(),
                "is_prokaryotic": _yes(r.get("is_prokaryotic")),
                "described": described,
                "described_source": (r.get("described_source") or "").strip(),
                "is_unrepresented": _yes(r.get("is_unrepresented")),
                "canonical_tree": (r.get("canonical_tree") or "").strip(),
                "has_tree": has_tree,
                "labels_raw": _int(r.get("labels_raw")),
                "labels_unique": labels,
                "c0_units": c0_units,
                "c0_unresolved_labels": _int(r.get("c0_unresolved_labels")),
                "time_scale": (r.get("time_scale") or "").strip(),
                "f_mol": f_mol,
                "f_low": _float(r.get("f_low"), f_mol),
                "f_high": _float(r.get("f_high"), f_mol),
                "f_basis": (r.get("f_basis") or "").strip(),
                "is_chronogram": chrono,
                "roottotip_cv": _float(r.get("roottotip_cv")),
                "interval_recoverable": recoverable,
                "interval_kind": (r.get("interval_kind") or "").strip(),
                "posterior_n": _int(r.get("posterior_n"), 0),
                "topology_assembly": (r.get("topology_assembly") or "").strip(),
                "C0": c0, "C1": c1, "C1d": c1d, "C1b": c1b,
            })
    if not rows:
        raise SystemExit(f"coverage basis is empty: {path}")
    return rows


def eukaryotic(rows):
    return [r for r in rows if not r["is_prokaryotic"]]


def prokaryotic(rows):
    return [r for r in rows if r["is_prokaryotic"]]


def canonical(rows):
    return [r for r in rows if r["has_tree"]]


def denominator(rows) -> int:
    """Described species summed over the given partitions.

    A column sum, never a literal -- which is precisely why the withdrawn
    23.7 % headline's implied denominator (2,084,409, present in no file)
    cannot be reconstructed here.
    """
    return sum(r["described"] for r in rows)


def rung_totals(rows) -> dict:
    """Numerator sums for the four rungs over the given partitions."""
    return {k: sum(r[k] for r in rows) for k in RUNGS}


def coverage_pct(rows, denom_rows=None) -> dict:
    """Rung percentages of the eukaryotic described-species denominator."""
    d = denominator(denom_rows if denom_rows is not None else rows)
    if d <= 0:
        raise SystemExit("denominator is zero -- refusing to divide")
    return {k: 100.0 * v / d for k, v in rung_totals(rows).items()}


def partition_pct(row) -> dict:
    """Per-partition rung percentages of that partition's described count."""
    d = row["described"]
    if d <= 0:
        return {k: 0.0 for k in RUNGS}
    return {k: 100.0 * row[k] / d for k in RUNGS}


# ── The generated headline file, used as a cross-check only ────────────────
def load_headlines(path: pathlib.Path | None = None) -> dict:
    with open(path or HEADLINE_JSON) as fh:
        return json.load(fh)


def check_against_headlines(headlines, *, numerators=None, percents=None,
                            denominator_euk=None, tol=1e-4, label="figure"):
    """Raise unless what we recomputed matches the deposited headline file.

    The two paths are genuinely independent -- we sum the basis CSV, the
    headline file was written by ``scripts/rebuild_atlas.py`` -- so agreement
    is evidence, not tautology.  Disagreement means the figure is about to
    draw a number the paper does not claim, and that must stop the build.
    """
    problems = []

    if denominator_euk is not None:
        want = headlines["denominators"]["eukaryotic"]["described"]
        if int(denominator_euk) != int(want):
            problems.append(
                f"eukaryotic denominator {denominator_euk:,} != headline {want:,}")

    for got_map, key, human in ((numerators, "numerator_sum", "numerator"),
                                (percents, "percent", "percent")):
        if not got_map:
            continue
        for rung, got in got_map.items():
            want = headlines["rungs"][rung][key]
            if not math.isclose(got, want, rel_tol=0, abs_tol=tol):
                problems.append(
                    f"{rung} {human} {got!r} != headline {want!r}")

    if problems:
        raise SystemExit(
            f"[{label}] recomputed values disagree with "
            f"{HEADLINE_JSON.name}:\n  - " + "\n  - ".join(problems) +
            "\nRefusing to render a figure the deposit does not support."
        )


def rung_label(headlines, rung) -> str:
    """'C0 NAMED -- 21.62 %' style label, definition text straight from JSON."""
    pct = headlines["rungs"][rung]["percent"]
    return f"{rung} {RUNG_SHORT[rung]} — {pct:.2f}%"


def rung_definition(headlines, rung) -> str:
    return headlines["rung_definitions"][rung]


# ── Per-tree provenance for the partition canonicals ───────────────────────
S5_PROVENANCE = INPUTS_DIR / "table_s5_tree_provenance.csv"
DATA_PROVENANCE = DATA_ROOT / "data_provenance.csv"

ARCHIVAL_CLASSES = ("undated", "point_dated", "hpd_recoverable")

_YEAR_RE = re.compile(r"\b(19[89]\d|20[0-4]\d)\b")


def load_canonical_provenance(basis, headlines):
    """One row per partition-canonical tree, with a resolved publication year.

    Reads the rebuilt ledger Paper/inputs/table_s5_tree_provenance.csv and
    keeps only rows naming a partition that the coverage basis says carries a
    canonical tree.  That data test -- rather than a hard-coded
    "(condamine_bulk)" string -- is what drops the bulk family-set row.

    Years are joined from Data/data_provenance.csv on DOI, falling back to a
    year embedded in the study label.  A canonical with no resolvable year
    raises: a tree that quietly disappears from a year axis is a worse defect
    than a build that stops.

    Returns ``(rows, skipped)`` where ``skipped`` lists the non-partition
    ledger rows, so a caller can disclose them rather than pretend they do
    not exist.
    """
    by_partition = {r["partition"]: r for r in canonical(basis)}

    years_by_doi = {}
    with open(DATA_PROVENANCE, newline="") as fh:
        for r in csv.DictReader(fh):
            doi = (r.get("doi") or "").strip().lower()
            year = _int(r.get("year"), 0)
            if doi and year:
                years_by_doi.setdefault(doi, year)

    rows, skipped = [], []
    with open(S5_PROVENANCE, newline="") as fh:
        for r in csv.DictReader(fh):
            name = (r.get("partition_name") or "").strip()
            study = (r.get("study") or "").strip()
            if name not in by_partition:
                skipped.append({"partition": name, "study": study})
                continue
            cls = (r.get("archival_class") or "").strip()
            if cls not in ARCHIVAL_CLASSES:
                raise SystemExit(
                    f"{name}: unknown archival_class {cls!r}; expected one of "
                    f"{ARCHIVAL_CLASSES}")
            doi = (r.get("doi") or "").strip().lower()
            year, source = years_by_doi.get(doi), "data_provenance/doi"
            if year is None:
                m = _YEAR_RE.search(study)
                if m:
                    year, source = int(m.group(1)), "study label"
            if year is None:
                raise SystemExit(
                    f"{name} ({study!r}): no publication year in "
                    f"{DATA_PROVENANCE.name} by DOI and none in the study "
                    f"label. Refusing to drop a canonical tree silently.")
            row = dict(by_partition[name])
            row.update({
                "study": study,
                "doi": doi,
                "year": year,
                "year_source": source,
                "archival_class": cls,
                "dating_code": (r.get("dating_code") or "").strip(),
            })
            rows.append(row)

    n_canon = headlines["trees"]["n_partition_canonical"]
    if len(rows) != n_canon:
        raise SystemExit(
            f"{len(rows)} canonical rows resolved from "
            f"{S5_PROVENANCE.name} but the deposit reports {n_canon} "
            f"partition-canonical trees")

    measured = sum(1 for r in rows if r["is_chronogram"])
    want_chrono = headlines["trees"]["n_partition_canonical_chronograms"]
    if measured != want_chrono:
        raise SystemExit(
            f"{measured} measured chronograms among the canonicals but the "
            f"deposit reports {want_chrono}")

    return rows, skipped


# ── Category ordering, taken from the data, not hard-coded ─────────────────
CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods",
                  "Other animals", "Microbes & protists"]


def category_rank(cat: str) -> int:
    try:
        return CATEGORY_ORDER.index(cat)
    except ValueError:
        return len(CATEGORY_ORDER)


def savefig(plt, outdir: pathlib.Path, stem: str, **kw):
    """Write PNG + PDF and report both paths."""
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    plt.savefig(png, dpi=300, bbox_inches="tight", **kw)
    plt.savefig(pdf, bbox_inches="tight", **kw)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    return png, pdf
