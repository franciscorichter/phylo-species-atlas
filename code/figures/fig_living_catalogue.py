"""The living catalogue: the monthly procedure that keeps the atlas current.

Drawn from the procedure as the code and the Methods state it (living/README.md,
living/watch.py, living/monthly.sh, living/curate.py, Methods "Continuous data
collection"):
  1 Watch    the previous month's deposits of the three data archives;
  2 Search   the whole history of the archives, group by group;
  both measured and screened by watch.py in one scheduled run (the timer in
  living/systemd/), which changes no count;
  3 Screen and curation   a tree becomes a candidate by the rules of D65;
  curate.py writes the measurable half of a dossier and a curator reads the
  source study; the verdict is recorded with its rule and a quotation;
  an adopted tree enters a rebuild (rebuild_atlas.py; verify_atlas.py re-derives
  every number without importing it), which writes the per-group table the next
  screen reads its reference trees from, and the data release;
  the author approves the release, which is then published as the article, the
  data deposit and the site; reproduce_paper.py in the release recomputes the
  article's numbers from the release's files.
Three kinds of mark: computed by code (filled ink), judged by reading the source
study (open ink), decided by the author (blue, afc.PERSON).  The figure prints
two quantities, the release tag (CITATION.cff) and its tree count
(headline_numbers.json); it runs in the paper's repository and in the release.
"""

import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

INK = "#1f2328"
PERSON = afc.PERSON   # a step decided by a person (atlas_figure_common)
SOFT = "#8b9098"
FAINT = "#a3a8b0"
LINE = "#e3e5e8"

# ── What the figure prints: the release tag and its tree count ─────────────
head = afc.load_headlines()
trees = head["trees"]
basis = afc.load_basis()
euk = afc.eukaryotic(basis)
afc.check_against_headlines(head, denominator_euk=afc.denominator(euk),
                            numerators=afc.rung_totals(afc.canonical(euk)),
                            label="figure_living_catalogue")
CFF = next(p for p in (afc.REPO_ROOT / "release" / "files" / "CITATION.cff",
                       afc.REPO_ROOT / "CITATION.cff") if p.is_file())
cff_version = re.search(r"^version: (\S+)", CFF.read_text(encoding="utf-8"), re.M).group(1)
release_tag = "v" + cff_version
META = afc.PAPER_ROOT / "shared" / "metadata_v2.tex"
if META.is_file():   # in the paper's repository the paper must cite this release
    cited = re.search(r"\\newcommand\{\\IdDataRelease\}\{([^}]+)\}",
                      META.read_text(encoding="utf-8")).group(1)
    if cited != release_tag:
        raise SystemExit(f"the paper cites {cited}; CITATION.cff is {cff_version}")


def fmt(n):
    return f"{n:,}"


# ── Geometry (millimetres at print size) ───────────────────────────────────
W, H = 180, 66
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")

Y = 35                           # the main line
LW = 1.0
NODE_R, GATE_R = 1.8, 2.5
X_ARCH = 13                      # the archives
X_IN, DY = 38, 8                 # watch (above) and search (below) nodes
X_SCREEN, X_CURATE, X_REBUILD = 66, 86, 110
X_REL, R_REL = 132, 10.5         # the release disc
X_GATE, X_PROD = 152, 162


def head(x, y, direction_deg, color, size=2.0):
    """A filled arrowhead with its tip at (x, y), pointing along direction."""
    t = math.radians(direction_deg)
    ux, uy = math.cos(t), math.sin(t)
    nx, ny = -uy, ux
    pts = [(x, y),
           (x - size * ux + 0.55 * size * nx, y - size * uy + 0.55 * size * ny),
           (x - size * ux - 0.55 * size * nx, y - size * uy - 0.55 * size * ny)]
    ax.add_patch(plt.Polygon(pts, closed=True, facecolor=color, edgecolor="none",
                             zorder=2))


def arrow(p0, p1, color=INK, gap0=0.0, gap1=0.0, lw=LW):
    """A straight arrow from p0 to p1, shortened by the gaps, head at p1."""
    (x0, y0), (x1, y1) = p0, p1
    d = math.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / d, (y1 - y0) / d
    a = (x0 + ux * gap0, y0 + uy * gap0)
    b = (x1 - ux * gap1, y1 - uy * gap1)
    ax.plot([a[0], b[0] - ux * 0.3], [a[1], b[1] - uy * 0.3], color=color, lw=lw,
            zorder=1, solid_capstyle="butt")
    head(b[0], b[1], math.degrees(math.atan2(uy, ux)), color)


def node(x, y, kind):
    """kind: 'code' (filled ink), 'read' (open ink), 'author' (blue)."""
    if kind == "code":
        ax.add_patch(Circle((x, y), NODE_R, facecolor=INK, edgecolor="white",
                            lw=0.8, zorder=3))
    elif kind == "read":
        ax.add_patch(Circle((x, y), NODE_R + 0.3, facecolor="white",
                            edgecolor=INK, lw=1.2, zorder=3))
    else:
        ax.add_patch(Circle((x, y), GATE_R, facecolor="white", edgecolor=PERSON,
                            lw=1.3, zorder=3))
        ax.add_patch(Circle((x, y), 1.1, facecolor=PERSON, edgecolor="none",
                            zorder=4))


def title(x, y, text, ha="center", va="bottom", color=INK):
    ax.text(x, y, text, ha=ha, va=va, fontsize=7.4, color=color, fontweight="bold")


def note(x, y, text, ha="center", va="top"):
    ax.text(x, y, text, ha=ha, va=va, fontsize=6.3, color=SOFT, linespacing=1.1)


# the three data archives (named in Methods)
ax.add_patch(plt.Rectangle((X_ARCH - 10, Y - 4.6), 20, 9.2, facecolor="white",
                           edgecolor=SOFT, lw=0.8, zorder=2))
ax.text(X_ARCH, Y, "data archives", ha="center", va="center", fontsize=6.8,
        color=INK, zorder=3)

# 1 Watch (the previous month) and 2 Search (the whole history), both measured
# and screened in one scheduled run
yw, ys = Y + DY, Y - DY
arrow((X_ARCH + 10, Y), (X_IN, yw), gap0=0.3, gap1=NODE_R + 0.7)
arrow((X_ARCH + 10, Y), (X_IN, ys), gap0=0.3, gap1=NODE_R + 0.7)
node(X_IN, yw, "code")
node(X_IN, ys, "code")
title(X_IN, yw + 3.6, "1  Watch")
title(X_IN, ys - 3.6, "2  Search", va="top")
arrow((X_IN, yw), (X_SCREEN, Y), gap0=NODE_R + 0.7, gap1=NODE_R + 0.7)
arrow((X_IN, ys), (X_SCREEN, Y), gap0=NODE_R + 0.7, gap1=NODE_R + 0.7)

# 3 Screen and curation: a computed screen, then a judged verdict
node(X_SCREEN, Y, "code")
node(X_CURATE, Y, "read")
arrow((X_SCREEN, Y), (X_CURATE, Y), gap0=NODE_R + 0.7, gap1=NODE_R + 1.0)
title((X_SCREEN + X_CURATE) / 2, Y + 4.8, "3  Screen and curation")

# the scheduled part: watch, search and screen run each month on a timer
from matplotlib.patches import FancyBboxPatch
xb0, xb1 = X_IN - 7.5, X_SCREEN + 4.6
yb0, yb1 = Y - DY - 7.7, Y + DY + 7.7
ax.add_patch(FancyBboxPatch((xb0, yb0), xb1 - xb0, yb1 - yb0,
                            boxstyle="round,pad=0,rounding_size=2.5",
                            facecolor="#f7f8f9", edgecolor=LINE, lw=0.9, zorder=0))
note((xb0 + xb1) / 2, yb1 + 0.8, "runs each month on a timer; changes no count",
     va="bottom")

# an adopted tree enters a rebuild, which writes the release
arrow((X_CURATE, Y), (X_REBUILD, Y), gap0=NODE_R + 1.0, gap1=NODE_R + 0.7)
note((X_CURATE + X_REBUILD) / 2, Y + 1.7, "adopted tree", va="bottom")
node(X_REBUILD, Y, "code")
title(X_REBUILD, Y + 4.8, "Rebuild")
arrow((X_REBUILD, Y), (X_REL - R_REL, Y), gap0=NODE_R + 0.7, gap1=0.3)

# the rebuild writes the per-group table from which the screen takes the
# reference trees of the next run
ax.add_patch(FancyArrowPatch((X_REBUILD - 0.8, Y - NODE_R - 0.8),
                             (X_SCREEN + 0.8, Y - NODE_R - 0.8),
                             arrowstyle="-|>", mutation_scale=8, color=SOFT, lw=0.9,
                             connectionstyle="arc3,rad=-0.40", zorder=1,
                             shrinkA=0, shrinkB=0))
note((X_SCREEN + X_REBUILD) / 2, Y - 15.2,
     "the per-group table: the reference trees of the next screen", va="top")

# the release
ax.add_patch(Circle((X_REL, Y), R_REL, facecolor=INK, edgecolor="none", zorder=2))
ax.text(X_REL, Y + 3.8, "data release", ha="center", va="center", fontsize=6.0,
        color=LINE, zorder=3)
ax.text(X_REL, Y - 0.2, release_tag, ha="center", va="center", fontsize=9.5,
        color="white", fontweight="bold", zorder=3)
ax.text(X_REL, Y - 4.4, f"{fmt(trees['n_shipped'])} trees", ha="center",
        va="center", fontsize=6.0, color="white", zorder=3)

# the author approves the release; then it is published
arrow((X_REL + R_REL, Y), (X_GATE, Y), gap0=0.3, gap1=GATE_R + 0.7)
node(X_GATE, Y, "author")
products = [(Y + 13, "Paper", 0.18), (Y, "Data deposit", 0.0), (Y - 13, "Site", -0.18)]
for py, text, rad in products:
    ax.add_patch(FancyArrowPatch((X_GATE + GATE_R + 0.4, Y), (X_PROD - 2.0, py),
                                 arrowstyle="-|>", mutation_scale=8,
                                 connectionstyle=f"arc3,rad={rad}", color=INK,
                                 lw=LW, zorder=2, shrinkA=0, shrinkB=0))
    ax.text(X_PROD, py, text, ha="left", va="center", fontsize=7.0, color=INK,
            fontweight="bold")

# key, bottom left
KY = 4.4
for k, (kind, txt) in enumerate([("code", "computed by code"),
                                 ("read", "judged by reading the source study"),
                                 ("author", "decided by the author")]):
    x = 6 + k * 46
    node(x, KY, kind)
    ax.text(x + 4, KY, txt, va="center", fontsize=6.5,
            color=PERSON if kind == "author" else INK)

afc.savefig(plt, afc.FIGURES_DIR, "figure_living_catalogue")
print(f"release {release_tag}, {trees['n_shipped']} trees  [headline_numbers.json, CITATION.cff]")
