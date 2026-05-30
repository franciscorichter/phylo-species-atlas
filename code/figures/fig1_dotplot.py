"""Figure 1: Per-partition phylogenetic coverage dotplot (print version).

Mirrors the live atlas page (franciscorichter.github.io/phylo-species-atlas).
For each of 62 partitions:
  ◇ filled marker  = species in shipped tree (tip count)
  ○ open circle    = described species
  ━━━━ horizontal interval = estimated true-diversity range (low–high)
Colour of the tip marker: blue = dated, amber = undated.
X-axis: log species count. Y-axis: partitions grouped by category band.

Data: data/processed_manuscript_inputs/fig1_coverage_bars.csv
Output: paper/figures/figure_coverage_dotplot.{png,pdf}
"""
import csv
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle

script_dir = pathlib.Path(__file__).resolve().parent
root = script_dir.parent.parent
csv_path = root / "inputs" / "fig1_coverage_bars.csv"
outdir = root / "figures"

mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

# ── Load ──
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

# ── Order: by category (Vertebrates → Microbes), then by described count desc ──
CAT_ORDER = ["Vertebrates", "Plants", "Arthropods", "Other animals", "Microbes & protists"]
CAT_RANK = {c: i for i, c in enumerate(CAT_ORDER)}
rows.sort(key=lambda r: (CAT_RANK.get(r["category"], 99),
                         -(r["estimated_total"] or r["described"] or 0)))

# ── Layout ──
n = len(rows)
fig_h = max(7.5, n * 0.18 + 1.5)
fig, ax = plt.subplots(figsize=(8.5, fig_h), constrained_layout=True)

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
DATED_COLOR = "#2a6fbf"
UNDATED_COLOR = "#c97a2a"
DESC_COLOR = "#222"
INTERVAL_COLOR = "#6b8cb3"
INTERVAL_HEURISTIC = "#d4b483"

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

    # Tip count — filled diamond (only if shipped tree)
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
    mlines.Line2D([], [], marker="D", color=DATED_COLOR, markersize=6, linewidth=0,
                  label="Shipped tree (dated)"),
    mlines.Line2D([], [], marker="D", color=UNDATED_COLOR, markersize=6, linewidth=0,
                  label="Shipped tree (undated)"),
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
