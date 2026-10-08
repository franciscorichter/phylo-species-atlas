"""Figure 4: crown ages of the reference chronograms and the divergence-time
interval each deposit carries, coloured by whether the interval counts for
interval coverage and, if not, why.

Input, written by scripts/rebuild_atlas.py:
  Paper/inputs/crown_ages.csv      one row per reference chronogram
  Paper/inputs/headline_numbers.json   cross-check of the counts
Output: Paper/figures/figure_crown_chronograms.{png,pdf}
"""
import csv
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

# colour, marker and legend label of each status; the marker separates the
# four classes without relying on hue (R2-39)
STATUS = OrderedDict([
    ("counted", ("#0072B2", "o", "interval usable")),
    ("kind", ("#E69F00", "s", "interval describes other variation")),
    ("fixed", ("#D55E00", "D", "excluded by a calibration rule")),
    ("none", ("#999999", "x", "no interval archived")),
])

rows = []
with open(afc.INPUTS_DIR / "crown_ages.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if not r["interval_kind"]:
            st = "none"
        elif r["counts_for_c1b"] == "yes":
            st = "counted"
        elif r["c1b_exclusion"] == "kind":
            st = "kind"
        else:
            st = "fixed"
        rows.append({
            "partition": r["partition"], "clade": r["crown_clade"],
            "age": float(r["crown_age_ma"]),
            "lo": float(r["interval_lo_ma"]) if r["interval_lo_ma"] else None,
            "hi": float(r["interval_hi_ma"]) if r["interval_hi_ma"] else None,
            "status": st,
        })

head = afc.load_headlines()
ca = head["crown_ages"]
if len(rows) != ca["n_absolute_chronograms"]:
    raise SystemExit("crown_ages.csv disagrees with headline_numbers.json")
if sum(r["status"] == "counted" for r in rows) != ca["n_with_counted_interval"]:
    raise SystemExit("counted intervals disagree with headline_numbers.json")

rows.sort(key=lambda r: r["age"])
fig, ax = plt.subplots(figsize=(5.2, 0.2 * len(rows) + 1.5),
                       constrained_layout=True)
for y, r in enumerate(rows):
    col, mk, _ = STATUS[r["status"]]
    if r["lo"] is not None:
        ax.plot([r["lo"], r["hi"]], [y, y], color=col, linewidth=2.2,
                solid_capstyle="butt", alpha=0.85)
    ax.plot(r["age"], y, marker=mk, color=col, markersize=4.6,
            markeredgecolor="white" if mk != "x" else col,
            markeredgewidth=0.5 if mk != "x" else 1.2, zorder=3,
            linestyle="none")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f"{r['partition']} ({r['clade'].split(' (')[0]})"
                    for r in rows], fontsize=6.3)
ax.set_xlabel("Crown age (Ma)")
ax.invert_xaxis()
ax.grid(axis="x", color="#e5e5e5", linewidth=0.5)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
handles = [mlines.Line2D([], [], color=c, marker=mk, linewidth=2.2
                         if mk != "x" else 0, markersize=4.6,
                         markeredgewidth=1.2 if mk == "x" else 0.5,
                         label=lab)
           for c, mk, lab in STATUS.values()]
# legend outside the axes, below the plot, so it covers no row
fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=6,
           frameon=False)

afc.savefig(plt, afc.FIGURES_DIR, "figure_crown_chronograms")
