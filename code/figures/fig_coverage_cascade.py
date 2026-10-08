"""Figure 1: described eukaryotic species at four levels of evidence, and
where coverage is lost between levels.

(a) Stacked bars: species named in a reference tree, placed with molecular
    data, in a chronogram in absolute time, and in a chronogram whose deposit
    carries a divergence-time interval that a calibration does not fix;
    stacked by group category.
(b) The loss between consecutive levels, split by reason.

Inputs, written by scripts/rebuild_atlas.py:
  Paper/inputs/fig1_coverage_cascade.csv   levels by category, losses by reason
  Paper/inputs/headline_numbers.json       cross-check of the totals
Output: Paper/figures/figure_coverage_cascade.{png,pdf}
"""
import csv
import math
from collections import OrderedDict

import matplotlib as mpl
import matplotlib.pyplot as plt

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

CAT_COLOR = {
    "Vertebrates": "#0072B2",
    "Plants": "#009E73",
    "Arthropods": "#D55E00",
    "Other animals": "#CC79A7",
    "Microbes & protists": "#E69F00",
}
LEVELS = OrderedDict([("C0", "Named"), ("C1", "Molecular"),
                      ("C1d", "Dated"), ("C1b", "Interval")])
REASONS = OrderedDict([
    ("TaxPlaced", ("Named to molecular",
                   "tips not in the molecular analysis")),
    ("NotChron", ("Molecular to dated", "tree is not a chronogram")),
    ("RelTime", ("Molecular to dated", "chronogram in relative time")),
    ("NoInterval", ("Dated to interval", "no interval archived")),
    ("Kind", ("Dated to interval", "interval describes other variation")),
    ("Calib", ("Dated to interval", "interval excluded by a calibration rule")),
])

level, loss = OrderedDict(), OrderedDict()
with open(afc.INPUTS_DIR / "fig1_coverage_cascade.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        v = float(r["species"])
        if r["kind"] == "level":
            level.setdefault(r["level"], OrderedDict())[r["group"]] = v
        else:
            loss[r["group"]] = v

head = afc.load_headlines()
denom = head["rungs"]["C0"]["denominator"]
for k in LEVELS:
    tot = sum(level[k].values())
    if tot != head["rungs"][k]["numerator_union"]:
        raise SystemExit(f"{k}: cascade total {tot} disagrees with "
                         f"headline_numbers.json")
if set(loss) != set(REASONS):
    raise SystemExit("cascade losses do not match the reasons drawn")

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1.0, 1.25]},
                             constrained_layout=True)
xs = range(len(LEVELS))
for i, k in enumerate(LEVELS):
    bottom = 0.0
    for cat in afc.CATEGORY_ORDER:
        v = level[k].get(cat, 0.0)
        ax.bar(i, v / 1e3, bottom=bottom / 1e3, width=0.66,
               color=CAT_COLOR[cat], edgecolor="white", linewidth=0.4,
               label=cat if i == 0 else None)
        bottom += v
    p = 100.0 * head["rungs"][k]["numerator_union"] / denom
    ax.text(i, bottom / 1e3, f"{head['rungs'][k]['numerator_union']:,}\n"
            f"({p:.1f}%)" if p >= 1 else
            f"{head['rungs'][k]['numerator_union']:,}\n({p:.2f}%)",
            ha="center", va="bottom", fontsize=6.5)
ax.set_xticks(list(xs))
ax.set_xticklabels(list(LEVELS.values()))
ax.set_ylabel("Described eukaryotic species (thousands)")
ymax = sum(level["C0"].values()) / 1e3
ax.set_ylim(0, ymax * 1.18)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(f"(a) Coverage of {denom:,} described species", loc="left",
             fontsize=7.5)
ax.legend(frameon=False, fontsize=6, loc="upper right")

ys = list(range(len(REASONS), 0, -1))
prev = None
for y, (k, (step, why)) in zip(ys, REASONS.items()):
    bx.barh(y, loss[k] / 1e3, height=0.62, color="#8c8c8c")
    bx.text(loss[k] / 1e3, y, f"  {loss[k]:,.0f}", va="center", fontsize=6.5)
bx.set_yticks(ys)
bx.set_yticklabels([f"{s}: {w}" for s, w in REASONS.values()], fontsize=6.5)
bx.set_xlabel("Species lost between levels (thousands)")
bx.set_xlim(0, max(loss.values()) / 1e3 * 1.3)
bx.spines[["top", "right"]].set_visible(False)
bx.set_title("(b) Where coverage is lost, by reason", loc="left",
             fontsize=7.5)

afc.savefig(plt, afc.FIGURES_DIR, "figure_coverage_cascade")
