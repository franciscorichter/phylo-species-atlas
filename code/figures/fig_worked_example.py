"""Figure 5: the worked example.  On each eligible reference chronogram,
(a) extinction fraction epsilon = mu / lambda on the molecular tree under the
sampling fraction that matches it (molecular tips / described) and under the
full tree's fraction (tips / described, decision D27); (b) net diversification r = lambda - mu on the full tree and on the
molecular tree, each at its own sampling fraction.  For a group whose source
archives a posterior set, else a set with dates drawn at random (column basis
of the table, decision D54), the point is the median and the bar the range
over the trees of that set (square markers); otherwise the point is the
summary-tree estimate and the bar its 95% profile-likelihood interval,
conditional on that one tree (round markers).

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
SET_MARK = {"posterior": "\u2020", "date_random": "\u2021"}
for y, r in enumerate(rows):
    f = lambda k: float(r[k])
    is_set = r["basis"] != "summary"
    mk = "s" if is_set else "o"
    for dy, k, c in ((off, "eps_molecular", C_MATCH),
                     (-off, "eps_fullrho", C_NAMED)):
        ax.plot([f(k + "_lo"), f(k + "_hi")], [y + dy] * 2, color=c,
                linewidth=1.6)
        ax.plot(f(k), y + dy, mk, color=c, markersize=3.6)
    if is_set and int(r["n_boundary"]):
        label = (f"\u03bc forced to zero in {r['n_boundary']} of {r['n_trees']} "
                 f"trees; LR median {f('lr'):.2f}")
    elif not is_set and r["classification"].startswith("boundary"):
        # D59: an outcome only on a binary summary tree
        tf = head["worked_example_S13"]["test_family"]
        why = {"polytomies": "polytomies",
               "node_ages_assigned": "node ages assigned"}
        mark = ("no test (%s)" % why[tf["not_tested"][r["partition"]]]
                if r["partition"] not in tf["tested"]
                else "rejected" if r["classification"].endswith("(rejected)")
                else "not rejected")
        label = f"\u03bc forced to zero; LR {f('lr'):.2f}, {mark} on this tree"
    else:
        label = ""
    if label:
        ax.annotate(label, (f("eps_fullrho_hi"), y - off), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=5.6,
                    color=C_NAMED)
    for dy, k, c in ((off, "r_full", C_FULL), (-off, "r_molecular", C_MOL)):
        bx.plot([f(k + "_lo"), f(k + "_hi")], [y + dy] * 2, color=c,
                linewidth=1.6)
        bx.plot(f(k), y + dy, mk, color=c, markersize=3.6)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r["partition"] + SET_MARK.get(r["basis"], "") for r in rows])
ax.set_xlabel(r"Extinction fraction $\varepsilon = \mu/\lambda$, pruned tree")
ax.set_xlim(left=-0.03)
bx.set_xlabel(r"Net diversification $r = \lambda - \mu$ (per Myr)")
for a in (ax, bx):
    a.grid(axis="x", color="#e5e5e5", linewidth=0.5)
    a.set_axisbelow(True)
    a.spines[["top", "right"]].set_visible(False)
ax.set_title("(a) Same tree, two sampling fractions", loc="left", fontsize=7.5)
bx.set_title("(b) Full against pruned tree", loc="left", fontsize=7.5)
ax.legend(handles=[
    mlines.Line2D([], [], color=C_MATCH, marker="o", label="sequenced tips / described"),
    mlines.Line2D([], [], color=C_NAMED, marker="o", label="full-tree tips / described")],
    loc="upper left", fontsize=6, frameon=False)
bx.legend(handles=[
    mlines.Line2D([], [], color=C_FULL, marker="o", label="full tree, tips / described"),
    mlines.Line2D([], [], color=C_MOL, marker="o",
                  label="pruned tree, sequenced tips / described")],
    loc="upper right", fontsize=6, frameon=False)

afc.savefig(plt, afc.FIGURES_DIR, "figure_worked_example")
