"""Figure 5: the worked example.  On each eligible reference chronogram,
(a) extinction fraction epsilon = mu / lambda on the molecular tree under the
sampling fraction that matches it (molecular tips / described) and under the
full tree's fraction (tips / described, decision D27); (b) net diversification r = lambda - mu on the full tree and on the
molecular tree, each at its own sampling fraction.  Bars are 95% profile-
likelihood intervals, conditional on the single tree fitted.

Input, written by scripts/rebuild_atlas.py:
  Paper/inputs/table_worked_example.csv
  Paper/inputs/headline_numbers.json   cross-check of the eligible groups
Output: Paper/figures/figure_worked_example.{png,pdf}
"""
import csv

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

rows = list(csv.DictReader(open(afc.INPUTS_DIR / "table_worked_example.csv",
                                newline="")))
head = afc.load_headlines()
if sorted(r["partition"] for r in rows) != sorted(
        head["worked_example_S13"]["eligible"]):
    raise SystemExit("table_worked_example.csv disagrees with the eligible set")
rows.sort(key=lambda r: float(r["eps_molecular"]))

C_MATCH, C_NAMED = "#0072B2", "#D55E00"
C_FULL, C_MOL = "#555555", "#009E73"
off = 0.17
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True,
                             constrained_layout=True)
for y, r in enumerate(rows):
    f = lambda k: float(r[k])
    for dy, k, c in ((off, "eps_molecular", C_MATCH),
                     (-off, "eps_fullrho", C_NAMED)):
        ax.plot([f(k + "_lo"), f(k + "_hi")], [y + dy] * 2, color=c,
                linewidth=1.6)
        ax.plot(f(k), y + dy, "o", color=c, markersize=3.6)
    if r["classification"].startswith("boundary"):
        mark = "rejected" if r["classification"].endswith("(rejected)") \
            else "not rejected"
        ax.annotate(f"\u03bc forced to zero; LR {f('lr'):.2f}, {mark}",
                    (f("eps_fullrho_hi"), y - off), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=5.6,
                    color=C_NAMED)
    for dy, k, c in ((off, "r_full", C_FULL), (-off, "r_molecular", C_MOL)):
        bx.plot([f(k + "_lo"), f(k + "_hi")], [y + dy] * 2, color=c,
                linewidth=1.6)
        bx.plot(f(k), y + dy, "o", color=c, markersize=3.6)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r["partition"] for r in rows])
ax.set_xlabel(r"Extinction fraction $\varepsilon = \mu/\lambda$, molecular tree")
ax.set_xlim(left=-0.03)
bx.set_xlabel(r"Net diversification $r = \lambda - \mu$ (per Myr)")
for a in (ax, bx):
    a.grid(axis="x", color="#e5e5e5", linewidth=0.5)
    a.set_axisbelow(True)
    a.spines[["top", "right"]].set_visible(False)
ax.set_title("(a) Same tree, two sampling fractions", loc="left", fontsize=7.5)
bx.set_title("(b) Full against molecular tree", loc="left", fontsize=7.5)
ax.legend(handles=[
    mlines.Line2D([], [], color=C_MATCH, marker="o", label="molecular tips / described"),
    mlines.Line2D([], [], color=C_NAMED, marker="o", label="full-tree tips / described")],
    loc="lower right", fontsize=6, frameon=False)
bx.legend(handles=[
    mlines.Line2D([], [], color=C_FULL, marker="o", label="full tree, tips / described"),
    mlines.Line2D([], [], color=C_MOL, marker="o",
                  label="molecular tree, molecular tips / described")],
    loc="upper right", fontsize=6, frameon=False)

afc.savefig(plt, afc.FIGURES_DIR, "figure_worked_example")
