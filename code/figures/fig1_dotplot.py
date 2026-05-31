"""Figure 1: Per-partition phylogenetic coverage dotplot (print version).

Mirrors the live atlas page (franciscorichter.github.io/phylo-species-atlas).
For each of 62 partitions:
  ◆ blue diamond   = atlas tree size (dated) — tip count of canonical dated tree
  ◆ green diamond  = molecular-only species (dated) — tips with direct molecular
                     sequence support (ntips × molecular_tip_fraction, from S5)
  ◆ amber diamond  = atlas tree size (undated)
  ○ open circle    = described species
  ━━━━ horizontal interval = estimated true-diversity range (low–high)
X-axis: log species count. Y-axis: partitions grouped by category band.

Inputs:
  - inputs/fig1_coverage_bars.csv          (62 rows, canonical plot source)
  - figure1-molecular-overlay-.../figure1_molecular_overlay.csv
        side-table keyed by partition_slug; joined on lowercase(group).
        Provides blue_diamond_x (canonical), green_diamond_x (molecular-only),
        amber_diamond_x (undated). For the 5 canonical eukaryote rows that
        appear as is_unrepresented=True in fig1_coverage_bars.csv, the overlay
        CSV is treated as the source of truth and the rows are rendered as
        dated with their blue+green diamonds.

Output: paper/figures/figure_coverage_dotplot.{png,pdf}
"""
import csv
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle

# ── Paths ──
# Canonical inputs live under the Paper tree. The two synced copies of this
# script (Code/.../code/figures and System/.../Paper/code/figures) both read
# the same absolute source so behaviour is identical regardless of which copy
# is invoked.
PAPER_ROOT = pathlib.Path("/Users/pancho/System/Research/phylo-atlas/Paper")
csv_path = PAPER_ROOT / "inputs" / "fig1_coverage_bars.csv"
overlay_csv_path = PAPER_ROOT / "figure1-molecular-overlay-2026-05-31" / "figure1_molecular_overlay.csv"
outdir = PAPER_ROOT / "figures"

mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

# ── Load main CSV ──
def parse_int(x):
    try: return int(float(x))
    except: return 0

rows = []
for r in csv.DictReader(open(csv_path)):
    rows.append({
        "group": r["group"],
        "category": r["category"],
        "tips": parse_int(r["tips"]) if r["tips"] else 0,
        "described": parse_int(r["described"]),
        "estimated_low": parse_int(r["estimated_low"]),
        "estimated_high": parse_int(r["estimated_high"]),
        "estimated_total": parse_int(r["estimated_total"]),
        "dated": r["dated"] == "True",
        "is_unrepresented": r["is_unrepresented"] == "True",
    })

# ── Load overlay CSV (side table) ──
# Schema: partition_slug, partition_name, blue_diamond_x, green_diamond_x,
#         amber_diamond_x, dated (yes/no)
def parse_overlay_int(x):
    if x is None or x == "" or x == "NA":
        return None
    try: return int(float(x))
    except: return None

overlay = {}  # key = lowercase partition_slug; value = dict
for r in csv.DictReader(open(overlay_csv_path)):
    slug = r["partition_slug"].strip().lower()
    overlay[slug] = {
        "blue":  parse_overlay_int(r.get("blue_diamond_x")),
        "green": parse_overlay_int(r.get("green_diamond_x")),
        "amber": parse_overlay_int(r.get("amber_diamond_x")),
        "dated": r.get("dated", "").strip().lower() == "yes",
    }

def overlay_key(group_name: str) -> str:
    """Match a fig1_coverage_bars 'group' value to an overlay partition_slug.

    Overlay slugs are lowercase with underscores ('seed_plants'); group names
    are title-case with spaces ('Seed plants'). Normalise both sides.
    """
    return group_name.strip().lower().replace(" ", "_")

# ── Reconcile: rows that are canonical-eukaryote dated in overlay but appear
# as is_unrepresented=True in fig1_coverage_bars.csv (Birds, Sharks, Seed
# plants, Squamates, Amphibians, ...). We trust the overlay: such a row should
# render as dated, with blue and green diamonds at the overlay values.
for r in rows:
    ov = overlay.get(overlay_key(r["group"]))
    if ov is None:
        continue
    if ov["dated"] and ov["blue"] is not None:
        # Promote to a represented, dated row.
        r["is_unrepresented"] = False
        r["dated"] = True
        # Preserve original tip count if present; otherwise fall back to overlay
        # blue value (so the diamond positioning matches even if main CSV was
        # empty).
        if not r["tips"]:
            r["tips"] = ov["blue"]
    # Stash overlay values on the row for the draw pass.
    r["overlay_blue"]  = ov["blue"]
    r["overlay_green"] = ov["green"]
    r["overlay_amber"] = ov["amber"]
    r["overlay_dated"] = ov["dated"]

# ── Order: by category (Vertebrates → Microbes), then by described count desc ──
CAT_ORDER = ["Vertebrates", "Plants", "Arthropods", "Other animals", "Microbes & protists"]
CAT_RANK = {c: i for i, c in enumerate(CAT_ORDER)}
rows.sort(key=lambda r: (CAT_RANK.get(r["category"], 99),
                         -(r["estimated_total"] or r["described"] or 0)))

# ── Layout ──
n = len(rows)
fig_h = max(7.5, n * 0.18 + 1.5)
# Widened from (8.5, fig_h) → (11.5, fig_h): ~35% wider so the new blue+green
# diamond pair has horizontal room without crowding the labels.
fig, ax = plt.subplots(figsize=(11.5, fig_h), constrained_layout=True)

y_positions = list(range(n, 0, -1))  # top row = first item

# ── Category background bands ──
CAT_COLORS = {
    "Vertebrates":          "#f4f7f9",  # very light blue
    "Plants":               "#f4f9f4",  # very light green
    "Arthropods":           "#fdf4ec",  # very light amber
    "Other animals":        "#f7f4f9",  # very light purple
    "Microbes & protists":  "#f9f4f4",  # very light red
}
# Compute band y-extents
bands = []
for cat in CAT_ORDER:
    idxs = [i for i, r in enumerate(rows) if r["category"] == cat]
    if not idxs: continue
    y_top = y_positions[min(idxs)] + 0.5
    y_bot = y_positions[max(idxs)] - 0.5
    bands.append((cat, y_bot, y_top))

# ── Draw bands ──
ax.set_xscale("log")
ax.set_xlim(0.5, 5e7)

for cat, y_bot, y_top in bands:
    ax.add_patch(Rectangle(
        (0.5, y_bot), 5e7, y_top - y_bot,
        facecolor=CAT_COLORS.get(cat, "#f8f8f8"),
        edgecolor="none", zorder=0))
    # Label on the right
    ax.text(5.5e7, (y_bot + y_top) / 2, cat,
            ha="left", va="center", fontsize=9, fontweight="bold",
            color="#666", rotation=0, clip_on=False)

# ── Per-row glyphs ──
DATED_COLOR    = "#1f77b4"   # blue  — atlas tree size (dated)
MOLECULAR_COLOR = "#2ca02c"  # green — molecular-only species (dated)
UNDATED_COLOR  = "#c97a2a"   # amber — atlas tree size (undated)
DESC_COLOR     = "#222"
INTERVAL_COLOR = "#6b8cb3"

for i, r in enumerate(rows):
    y = y_positions[i]

    # Interval line (low → high)
    if r["estimated_low"] and r["estimated_high"]:
        ax.plot([r["estimated_low"], r["estimated_high"]], [y, y],
                color=INTERVAL_COLOR, linewidth=1.2, alpha=0.55, zorder=1)

    # Described species — open circle
    if r["described"]:
        ax.plot(r["described"], y, marker="o", markersize=5,
                markerfacecolor="white", markeredgecolor=DESC_COLOR,
                markeredgewidth=1.0, zorder=3)

    ov_blue  = r.get("overlay_blue")
    ov_green = r.get("overlay_green")
    ov_amber = r.get("overlay_amber")

    if ov_blue is not None or ov_green is not None or ov_amber is not None:
        # Overlay-driven row — use overlay values verbatim.
        if ov_blue is not None:
            ax.plot(ov_blue, y, marker="D", markersize=6.5,
                    markerfacecolor=DATED_COLOR, markeredgecolor=DATED_COLOR,
                    zorder=4)
        if ov_green is not None:
            # Slightly smaller + thin black edge so it remains visible when it
            # coincides with the blue diamond (molecular_tip_fraction = 1.00).
            ax.plot(ov_green, y, marker="D", markersize=4.5,
                    markerfacecolor=MOLECULAR_COLOR, markeredgecolor="black",
                    markeredgewidth=0.4, zorder=5)
        if ov_amber is not None:
            ax.plot(ov_amber, y, marker="D", markersize=5.5,
                    markerfacecolor=UNDATED_COLOR, markeredgecolor=UNDATED_COLOR,
                    zorder=4)
    else:
        # Fallback: rows not in the overlay table — preserve original behaviour.
        if r["tips"] and not r["is_unrepresented"]:
            color = DATED_COLOR if r["dated"] else UNDATED_COLOR
            ax.plot(r["tips"], y, marker="D", markersize=5.5,
                    markerfacecolor=color, markeredgecolor=color,
                    zorder=4)

# ── Y-axis: partition labels ──
ax.set_yticks(y_positions)
ax.set_yticklabels([r["group"] for r in rows], fontsize=7)
ax.set_ylim(0, n + 1)
ax.tick_params(axis="y", length=0, pad=2)

# ── X-axis ──
ax.set_xlabel("Species count (log scale)", fontsize=9)
ax.set_xticks([1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7])
ax.set_xticklabels(["1", "10", "100", "1k", "10k", "100k", "1M", "10M"])

# ── Legend (top of plot) ──
import matplotlib.lines as mlines
legend_elements = [
    mlines.Line2D([], [], marker="D", color=DATED_COLOR, markersize=7, linewidth=0,
                  label="Blue diamond = atlas tree size (dated)"),
    mlines.Line2D([], [], marker="D", color=MOLECULAR_COLOR, markersize=5,
                  markeredgecolor="black", markeredgewidth=0.4, linewidth=0,
                  label="Green diamond = molecular-only species (dated)"),
    mlines.Line2D([], [], marker="D", color=UNDATED_COLOR, markersize=6, linewidth=0,
                  label="Amber diamond = atlas tree size (undated)"),
    mlines.Line2D([], [], marker="o", color=DESC_COLOR, markerfacecolor="white",
                  markersize=6, linewidth=0, label="Described species"),
    mlines.Line2D([], [], color=INTERVAL_COLOR, linewidth=2,
                  label="Estimated true-diversity interval"),
]
ax.legend(handles=legend_elements, loc="lower right", frameon=False,
          fontsize=8, ncol=1, bbox_to_anchor=(0.98, 0.02))

ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5, alpha=0.7, zorder=0)
ax.set_axisbelow(True)

# ── Save ──
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / "figure_coverage_dotplot.png", dpi=300, bbox_inches="tight")
plt.savefig(outdir / "figure_coverage_dotplot.pdf", bbox_inches="tight")
print(f"Wrote {outdir / 'figure_coverage_dotplot.png'} ({n} partitions)")
print(f"Wrote {outdir / 'figure_coverage_dotplot.pdf'}")
