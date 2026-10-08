"""Figure 2: named and molecular coverage, per eukaryotic group.

A diamond in panel (b) marks a group whose reference chronogram's deposit
carries a counted divergence-time interval (interval coverage > 0).

Two panels over the same eukaryotic partitions, sorted by C0:

  (a) C0  percent of the partition's described species that resolve, as tips
          of its canonical tree, to distinct Catalogue of Life accepted names
  (b) C1  the same quantity weighted by the measured molecular fraction,
          in tip-equivalents

Bars are solid where the canonical tree is a measured chronogram in absolute
time (so the partition contributes at C1d) and hatched where it is not.
Partitions with no canonical tree are structural zeros.  They are not drawn
(finding R1-34: 36 empty rows took most of the figure); the caption gives
their number and share of described species, and the rows of the basis file
are still cross-checked below.
The horizontal axis is logarithmic because per-partition coverage spans more
than three orders of magnitude.

Inputs, written by scripts/rebuild_atlas.py:
  Paper/inputs/fig3_c0c1_bars.csv          the per-partition projection
  Data/coverage_basis.csv                   cross-check
  Paper/inputs/headline_numbers.json        cross-check of the totals

Output: Paper/figures/figure_c0c1_bars.{png,pdf}
"""
import csv
import math

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

LEVEL_NAME = {"C0": "Named", "C1": "Molecular"}
SRC = afc.INPUTS_DIR / "fig3_c0c1_bars.csv"
rows = []
with open(SRC, newline="") as fh:
    for r in csv.DictReader(fh):
        rows.append({
            "partition": r["partition"],
            "category": r["category"],
            "described": afc._int(r["described"]),
            "c0_units": afc._int(r["c0_units"]),
            "C0": afc._float(r["coverage_C0_pct"], 0.0),
            "C1": afc._float(r["coverage_C1_pct"], 0.0),
            "zero": afc._yes(r["is_structural_zero"]),
            "chrono": afc._yes(r["is_chronogram"]),
        })

# ── Cross-checks ─────────────────────────────────────────────────────────────
basis = afc.load_basis()
head = afc.load_headlines()
euk = afc.eukaryotic(basis)
by_part = {b["partition"]: b for b in euk}
if sorted(by_part) != sorted(r["partition"] for r in rows):
    raise SystemExit("fig3_c0c1_bars.csv must list exactly the eukaryotic "
                     "partitions of Data/coverage_basis.csv")
for r in rows:
    b = by_part[r["partition"]]
    want0 = 100.0 * b["C0"] / b["described"] if b["described"] else 0.0
    want1 = 100.0 * b["C1"] / b["described"] if b["described"] else 0.0
    if not (math.isclose(r["C0"], want0, abs_tol=1e-3)
            and math.isclose(r["C1"], want1, abs_tol=1e-3)):
        raise SystemExit(f"{r['partition']}: bars input disagrees with the "
                         f"basis file")
    if r["zero"] == b["has_tree"]:
        raise SystemExit(f"{r['partition']}: structural-zero flag disagrees "
                         f"with the basis file")
    r["time_scale"] = b["time_scale"]
    r["interval"] = b["C1b"] > 0
afc.check_against_headlines(head, numerators=afc.rung_totals(euk),
                            denominator_euk=afc.denominator(euk),
                            label="fig_c0c1_bars")

# ── Order: groups with a reference tree, by C0 descending ───────────────────
n_zero = sum(r["zero"] for r in rows)
desc_zero = sum(r["described"] for r in rows if r["zero"])
rows = [r for r in rows if not r["zero"]]
rows.sort(key=lambda r: (-r["C0"], -r["described"]))
n = len(rows)
ys = list(range(n, 0, -1))

positive = [v for r in rows for v in (r["C0"], r["C1"]) if v > 0]
xmin = 10 ** math.floor(math.log10(min(positive)))
xmax = 10 ** math.ceil(math.log10(max(positive)))

fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.16 * n + 1.3), sharey=True,
                         constrained_layout=True)
denom = afc.denominator(euk)

for ax, rung, tag in ((axes[0], "C0", "a"), (axes[1], "C1", "b")):
    col = afc.RUNG_COLOR[rung]
    ax.set_xscale("log")
    ax.set_xlim(xmin, xmax * 1.6)
    for r, y in zip(rows, ys):
        if r["zero"]:
            ax.barh(y, xmax - xmin, left=xmin, height=0.8, color="#eeeeee",
                    edgecolor="#bbbbbb", hatch="////", linewidth=0.0,
                    zorder=1)
            continue
        v = r[rung]
        if v <= 0:
            continue
        ax.barh(y, v - xmin, left=xmin, height=0.8,
                color=col if r["chrono"] else "white",
                edgecolor=col, hatch=None if r["chrono"] else "....",
                linewidth=0.8, zorder=2)
        if rung == "C1" and r["interval"]:
            ax.plot(v * 1.3, y, marker="D", color="#222222", markersize=3.2,
                    zorder=3)
    ax.set_title(f"({tag}) {LEVEL_NAME[rung]} coverage, "
                 f"{head['rungs'][rung]['percent_union']:.1f}% overall",
                 loc="left", fontsize=7.5)
    # plain decimal tick labels: log-axis mathtext fails under Helvetica
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel("% of the group's described species (log scale)")
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)

axes[0].set_yticks(ys)
axes[0].set_yticklabels([r["partition"] for r in rows], fontsize=6.5)
axes[0].set_ylim(0.3, n + 0.7)

handles = [
    mpatches.Patch(facecolor=afc.RUNG_COLOR["C0"], edgecolor=afc.RUNG_COLOR["C0"],
                   label="named; reference tree is a chronogram in absolute time"),
    mpatches.Patch(facecolor=afc.RUNG_COLOR["C1"], edgecolor=afc.RUNG_COLOR["C1"],
                   label="molecular; reference tree is a chronogram in absolute time"),
    mpatches.Patch(facecolor="white", edgecolor="#555555", hatch="....",
                   label="reference tree is not a chronogram in absolute time"),
    plt.Line2D([], [], marker="D", color="#222222", linestyle="none",
               markersize=3.2,
               label="source archive carries a usable divergence-time interval"),
]
fig.legend(handles=handles, loc="outside lower center", ncol=2, fontsize=6.5,
           frameon=False)

afc.savefig(plt, afc.FIGURES_DIR, "figure_c0c1_bars")
print(f"rows: {n}; structural zeros: {n_zero} ({desc_zero:,} described)")
