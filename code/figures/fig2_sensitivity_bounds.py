"""Figure 2: Per-partition sensitivity bounds (C0 / C1 / C2).

For each partition we plot three coverage measures on a single horizontal
track so the reader can see how each partition's headline number drifts as
we tighten the inclusion criterion:

  C0 (permissive)  = filled blue diamond
                     all retained canonical tips counted against the
                     described-species denominator.
  C1 (strict)      = filled orange diamond
                     molecular-only tips on TACT-imputed trees, full count
                     elsewhere (i.e. demote heavily-imputed partitions).
  C2 (strictest)   = open red circle (or filled red X if the partition is
                     excluded entirely under C2 — e.g. full-OTL-pipeline
                     trees such as birds and insects).

A thin light-gray line connects C0 -> C1 -> C2 so the drift is visible at
a glance. Rows are sorted by drift_magnitude descending (most-affected
partitions on top) so the eye lands on the largest gaps first.

Inputs:
  figure2-3-rework-2026-05-31/c0c1c2_per_partition.csv

Outputs:
  figures/figure_sensitivity_bounds.{png,pdf}

Headline-anchor verticals:
  C0 global = 23.7%   (dashed blue)
  C1 global = 10.0%   (dashed orange)
  C2 global = 6.6%    (dashed red)
Per-partition values differ substantially from these globals (which are
weighted by described-species across the whole atlas); the verticals are
visual anchors, not per-partition baselines.
"""
import csv
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle

# ── Paths ──
PAPER_ROOT = pathlib.Path("/Users/pancho/System/Research/phylo-atlas/Paper")
csv_path = PAPER_ROOT / "figure2-3-rework-2026-05-31" / "c0c1c2_per_partition.csv"
outdir = PAPER_ROOT / "figures"

# ── Style ──
mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = 8
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

# ── Colours ──
C0_COLOR = "#1f77b4"   # blue   — permissive
C1_COLOR = "#ff7f0e"   # orange — strict molecular-only on TACT-imputed
C2_COLOR = "#d62728"   # red    — strictest, non-imputed + non-OTL
DRIFT_LINE_COLOR = "#bbbbbb"
EXCLUDED_C2_THRESHOLD = 1e-6   # treat coverage_c2_pct ≈ 0 as "excluded"

# ── Headline anchors (atlas-wide) ──
HEADLINE_C0 = 23.7
HEADLINE_C1 = 10.0
HEADLINE_C2 = 6.6

# ── Category display ──
CAT_DISPLAY = {
    "vertebrates":         "Vertebrates",
    "plants":              "Plants",
    "arthropods":          "Arthropods",
    "other invertebrates": "Other invertebrates",
    "protists":            "Protists & microbes",
}
CAT_ORDER = ["vertebrates", "plants", "arthropods",
             "other invertebrates", "protists"]
CAT_COLORS = {
    "vertebrates":         "#f4f7f9",
    "plants":              "#f4f9f4",
    "arthropods":          "#fdf4ec",
    "other invertebrates": "#f7f4f9",
    "protists":            "#f9f4f4",
}

# ── Load ──
def parse_float(x):
    try: return float(x)
    except: return 0.0

rows = []
with open(csv_path, newline="") as fh:
    for r in csv.DictReader(fh):
        if not r.get("partition_slug"):
            continue
        rows.append({
            "slug":     r["partition_slug"].strip(),
            "name":     r["partition_name"].strip(),
            "category": r["category"].strip().lower(),
            "c0":       parse_float(r["coverage_c0_pct"]),
            "c1":       parse_float(r["coverage_c1_pct"]),
            "c2":       parse_float(r["coverage_c2_pct"]),
            "drift":    parse_float(r["drift_magnitude"]),
        })

if not rows:
    raise SystemExit(f"No rows parsed from {csv_path}")

# ── Sort: by category in canonical order, then by drift_magnitude desc within ──
# This keeps the category bands intact but pushes the most-affected partitions
# to the top of each band — the eye still lands on the biggest drifts first.
cat_rank = {c: i for i, c in enumerate(CAT_ORDER)}
rows.sort(key=lambda r: (cat_rank.get(r["category"], 99),
                         -r["drift"], -r["c0"]))

n = len(rows)

# ── Identify top-3 drift partitions for annotation ──
top_drift_slugs = {r["slug"] for r in
                   sorted(rows, key=lambda r: -r["drift"])[:3]}

# ── Layout ──
fig_h = max(7.0, n * 0.34 + 1.4)
fig, ax = plt.subplots(figsize=(11.5, fig_h), constrained_layout=True)

y_positions = list(range(n, 0, -1))   # top row = first item

# ── Category background bands ──
bands = []
for cat in CAT_ORDER:
    idxs = [i for i, r in enumerate(rows) if r["category"] == cat]
    if not idxs:
        continue
    y_top = y_positions[min(idxs)] + 0.5
    y_bot = y_positions[max(idxs)] - 0.5
    bands.append((cat, y_bot, y_top))

ax.set_xlim(-2, 108)

for cat, y_bot, y_top in bands:
    ax.add_patch(Rectangle(
        (-2, y_bot), 110, y_top - y_bot,
        facecolor=CAT_COLORS.get(cat, "#f8f8f8"),
        edgecolor="none", zorder=0))
    ax.text(101.5, (y_bot + y_top) / 2, CAT_DISPLAY.get(cat, cat),
            ha="left", va="center", fontsize=9, fontweight="bold",
            color="#666", rotation=0, clip_on=False)

# ── Headline-anchor vertical reference lines ──
for x, color, label in [
    (HEADLINE_C0, C0_COLOR, f"C0 global = {HEADLINE_C0}%"),
    (HEADLINE_C1, C1_COLOR, f"C1 global = {HEADLINE_C1}%"),
    (HEADLINE_C2, C2_COLOR, f"C2 global = {HEADLINE_C2}%"),
]:
    ax.axvline(x, color=color, linestyle="--", linewidth=0.8,
               alpha=0.45, zorder=1)
    ax.text(x, n + 0.7, label, color=color, fontsize=7,
            ha="center", va="bottom", alpha=0.85, clip_on=False)

# ── Per-row glyphs ──
for i, r in enumerate(rows):
    y = y_positions[i]
    c0, c1, c2 = r["c0"], r["c1"], r["c2"]

    # Drift connector: only meaningful when the three points differ.
    xs = [c0, c1, c2]
    if max(xs) - min(xs) > 1e-6:
        ax.plot(xs, [y, y, y], color=DRIFT_LINE_COLOR,
                linewidth=1.0, alpha=0.85, zorder=2)

    # C0 — filled blue diamond
    ax.plot(c0, y, marker="D", markersize=7.0,
            markerfacecolor=C0_COLOR, markeredgecolor=C0_COLOR,
            zorder=4)

    # C1 — filled orange diamond
    ax.plot(c1, y, marker="D", markersize=6.0,
            markerfacecolor=C1_COLOR, markeredgecolor=C1_COLOR,
            markeredgewidth=0.4, zorder=5)

    # C2 — open red circle, OR filled red X if excluded under C2
    if c2 <= EXCLUDED_C2_THRESHOLD:
        ax.plot(0, y, marker="X", markersize=8.5,
                markerfacecolor=C2_COLOR, markeredgecolor=C2_COLOR,
                zorder=6)
    else:
        ax.plot(c2, y, marker="o", markersize=7.0,
                markerfacecolor="white", markeredgecolor=C2_COLOR,
                markeredgewidth=1.3, zorder=6)

    # Annotate top-3 drift partitions with their drift magnitude.
    if r["slug"] in top_drift_slugs and r["drift"] > 1.0:
        # Place the annotation to the right of the rightmost glyph.
        right_x = max(c0, c1, c2 if c2 > EXCLUDED_C2_THRESHOLD else 0)
        ax.annotate(f"Δ = {r['drift']:.1f} pp",
                    xy=(right_x, y),
                    xytext=(6, 0), textcoords="offset points",
                    fontsize=7, color="#444", va="center",
                    fontweight="bold")

# ── Y-axis: partition labels ──
ax.set_yticks(y_positions)
ax.set_yticklabels([r["name"] for r in rows], fontsize=8)
ax.set_ylim(0, n + 1)
ax.tick_params(axis="y", length=0, pad=2)

# ── X-axis ──
ax.set_xlabel("Coverage of described species (%)", fontsize=9)
ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax.grid(axis="x", which="major", color="#dddddd", linewidth=0.5,
        alpha=0.7, zorder=0)
ax.set_axisbelow(True)

# ── Legend ──
legend_elements = [
    mlines.Line2D([], [], marker="D", color=C0_COLOR, markersize=7,
                  linewidth=0,
                  label="C0 (permissive): all retained canonical tips"),
    mlines.Line2D([], [], marker="D", color=C1_COLOR, markersize=6,
                  linewidth=0,
                  label="C1 (strict): molecular-only on TACT-imputed trees"),
    mlines.Line2D([], [], marker="o", color=C2_COLOR,
                  markerfacecolor="white", markeredgewidth=1.3,
                  markersize=7, linewidth=0,
                  label="C2 (strictest): non-imputed, non-OTL-pipeline"),
    mlines.Line2D([], [], marker="X", color=C2_COLOR, markersize=8,
                  linewidth=0,
                  label="C2 = excluded entirely (full-OTL-pipeline tree)"),
    mlines.Line2D([], [], color=DRIFT_LINE_COLOR, linewidth=1.2,
                  label="Drift C0 -> C1 -> C2"),
]
ax.legend(handles=legend_elements, loc="lower right", frameon=False,
          fontsize=7.5, ncol=1, bbox_to_anchor=(0.985, 0.015))

# ── Save ──
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / "figure_sensitivity_bounds.png",
            dpi=300, bbox_inches="tight")
plt.savefig(outdir / "figure_sensitivity_bounds.pdf",
            bbox_inches="tight")
print(f"Wrote {outdir / 'figure_sensitivity_bounds.png'} "
      f"({n} partitions, top-drift = "
      f"{sorted(rows, key=lambda r: -r['drift'])[0]['name']})")
print(f"Wrote {outdir / 'figure_sensitivity_bounds.pdf'}")
