"""Figure: the two measured axes, shown together and not merged.

One point per eukaryotic partition-canonical tree:

  x  measured molecular fraction f_mol (with the band f_low..f_high where the
     deposit records one)
  y  C0 coverage, percent of that partition's own described-species count
     (log scale: the values span three orders of magnitude)
  marker area   proportional to described species
  marker shape  nominal topology_assembly annotation (never a criterion)
  fill          filled when the shipped canonical is a measured chronogram in
                absolute time (contributes at C1d); hollow when it is not

Inputs, both written by scripts/rebuild_atlas.py:
  Paper/inputs/fig2_coverage_scatter.csv   the per-partition projection
  Paper/inputs/headline_numbers.json       cross-check of the rung totals

No value is written down here.  Output:
  Paper/figures/figure_coverage_scatter.{png,pdf}
"""
import csv
import math

import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=8)

SRC = afc.INPUTS_DIR / "fig2_coverage_scatter.csv"

rows = []
with open(SRC, newline="") as fh:
    for r in csv.DictReader(fh):
        rows.append({
            "partition": r["partition"],
            "category": r["category"],
            "described": afc._int(r["described"]),
            "c0_units": afc._int(r["c0_units"]),
            "f_mol": afc._float(r["f_mol"]),
            "f_low": afc._float(r["f_low"]),
            "f_high": afc._float(r["f_high"]),
            "c0_pct": afc._float(r["coverage_C0_pct"]),
            "c1_pct": afc._float(r["coverage_C1_pct"]),
            "assembly": r["topology_assembly"],
            "chrono": afc._yes(r["is_chronogram"]),
            "label": r["label"],
        })

# ── Cross-check against the basis and the headline file ────────────────────
basis = afc.load_basis()
head = afc.load_headlines()
euk_canon = [b for b in afc.canonical(afc.eukaryotic(basis))]
if sorted(r["partition"] for r in rows) != sorted(b["partition"] for b in euk_canon):
    raise SystemExit("fig2_coverage_scatter.csv does not list exactly the "
                     "eukaryotic partition-canonical trees of the basis file")
by_part = {b["partition"]: b for b in euk_canon}
for r in rows:
    b = by_part[r["partition"]]
    if r["c0_units"] != b["c0_units"] or not math.isclose(
            r["f_mol"], b["f_mol"], abs_tol=1e-6):
        raise SystemExit(f"{r['partition']}: scatter input disagrees with "
                         f"Data/coverage_basis.csv")
afc.check_against_headlines(
    head,
    numerators={k: v for k, v in afc.rung_totals(afc.eukaryotic(basis)).items()
                if k in ("C0", "C1")},
    denominator_euk=afc.denominator(afc.eukaryotic(basis)),
    label="fig_coverage_scatter",
)

# ── Encoding ────────────────────────────────────────────────────────────────
SHAPE = {"primary_analysis": "o", "compiled": "s", "automated_synthesis": "^"}
SHAPE_LABEL = {"primary_analysis": "primary analysis",
               "compiled": "compiled from public sequences",
               "automated_synthesis": "automated synthesis"}
unknown = {r["assembly"] for r in rows} - set(SHAPE)
if unknown:
    raise SystemExit(f"unknown topology_assembly value(s): {sorted(unknown)}")

CAT_COLOR = {
    "Vertebrates": "#0072B2",
    "Plants": "#009E73",
    "Arthropods": "#D55E00",
    "Other animals": "#CC79A7",
    "Microbes & protists": "#E69F00",
}
max_desc = max(r["described"] for r in rows)
AREA_MAX = 600.0  # points^2 of the largest partition (a drawing width)


def area(d):
    return AREA_MAX * d / max_desc + 12.0


fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
ax.set_yscale("log")
ymin = min(r["c0_pct"] for r in rows)
ax.set_ylim(ymin / 2.5, 250)
ax.set_xlim(-0.03, 1.12)

for r in sorted(rows, key=lambda r: -r["described"]):
    col = CAT_COLOR.get(r["category"], "#555555")
    if r["f_low"] is not None and r["f_high"] is not None and \
            r["f_high"] - r["f_low"] > 1e-9:
        ax.plot([r["f_low"], r["f_high"]], [r["c0_pct"]] * 2, color=col,
                linewidth=2.0, alpha=0.6, zorder=2)
    ax.scatter(r["f_mol"], r["c0_pct"], s=area(r["described"]),
               marker=SHAPE[r["assembly"]],
               facecolor=col if r["chrono"] else "white",
               edgecolor=col, linewidth=1.1, alpha=0.85, zorder=3)

# Labels.  Points at f_mol = 1 alternate left/right in order of C0; every
# other label sits to the right.  Labels on the same side whose points are
# close in f_mol are then pushed apart vertically (in log10 C0) until no two
# are closer than LABEL_GAP decades, so names never print over one another.
LABEL_GAP = 0.085   # decades of the y axis between label centres
X_NEAR = 0.08       # f_mol distance within which labels can collide
at_one = sorted([r for r in rows if r["f_mol"] > 0.95],
                key=lambda r: r["c0_pct"])
for i, r in enumerate(at_one):
    r["right"] = (i % 2 == 1)
for r in rows:
    r.setdefault("right", True)
    r["ly"] = math.log10(r["c0_pct"])
for side in (True, False):
    lab = sorted([r for r in rows if r["right"] == side],
                 key=lambda r: r["ly"])
    for _ in range(200):
        moved = False
        for i in range(len(lab)):
            for j in range(i + 1, len(lab)):
                p_, q_ = lab[i], lab[j]
                if abs(p_["f_mol"] - q_["f_mol"]) > X_NEAR:
                    continue
                gap = q_["ly"] - p_["ly"]
                if abs(gap) < LABEL_GAP:
                    push = (LABEL_GAP - abs(gap)) / 2 + 1e-4
                    if gap >= 0:
                        p_["ly"] -= push
                        q_["ly"] += push
                    else:
                        p_["ly"] += push
                        q_["ly"] -= push
                    moved = True
        if not moved:
            break
for r in rows:
    right = r["right"]
    # clear the largest marker the label could overlap: any point within
    # X_NEAR in f_mol and LABEL_GAP decades of the label (R2-21)
    pad = max(math.sqrt(area(q["described"])) / 2 for q in rows
              if abs(q["f_mol"] - r["f_mol"]) <= X_NEAR
              and abs(math.log10(q["c0_pct"]) - r["ly"]) <= 1.5 * LABEL_GAP
              ) + 3   # marker radius + gap, pt
    ax.annotate(r["label"], (r["f_mol"], 10 ** r["ly"]),
                xytext=(pad if right else -pad, 0), textcoords="offset points",
                ha="left" if right else "right", va="center",
                fontsize=6.5, color="#333333", zorder=4)

ax.set_xlabel("Molecular fraction (molecular tips / tips of the reference tree)")
ax.set_ylabel("Named coverage (% of the group's described species, log scale)")
ax.set_yticks([0.1, 1, 10, 100])
ax.set_yticklabels(["0.1", "1", "10", "100"])
ax.grid(color="#e5e5e5", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

handles = []
for cat, col in CAT_COLOR.items():
    if any(r["category"] == cat for r in rows):
        handles.append(mlines.Line2D([], [], marker="o", linewidth=0,
                                     color=col, markersize=6, label=cat))
for key, mk in SHAPE.items():
    if any(r["assembly"] == key for r in rows):
        handles.append(mlines.Line2D([], [], marker=mk, linewidth=0,
                                     color="#555555", markersize=6,
                                     label=SHAPE_LABEL[key]))
handles.append(mlines.Line2D([], [], marker="o", linewidth=0,
                             color="#555555", markersize=6,
                             label="chronogram in absolute time (filled)"))
handles.append(mlines.Line2D([], [], marker="o", linewidth=0,
                             markerfacecolor="white", color="#555555",
                             markersize=6, label="not a chronogram in absolute time (hollow)"))
ax.legend(handles=handles, loc="lower left", fontsize=6.5, frameon=False,
          ncol=2)

afc.savefig(plt, afc.FIGURES_DIR, "figure_coverage_scatter")
print(f"points: {len(rows)}; measured chronograms: "
      f"{sum(r['chrono'] for r in rows)}")
