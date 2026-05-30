"""Figure 4 (replaces global LTT): Atlas growth through time.

Shows how the atlas would have looked at each historical year, given that
each tree could only have been included after its publication year:
  - cumulative number of trees
  - cumulative number of partitions covered (at least one canonical tree)
  - cumulative unique species in canonical partition trees

Data: data/processed_manuscript_inputs/figN_atlas_growth.csv
Output: paper/figures/figure_atlas_growth.{png,pdf}
"""
import csv
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl

script_dir = pathlib.Path(__file__).resolve().parent
root = script_dir.parent.parent
data_csv = root / "inputs" / "figN_atlas_growth.csv"
outdir = root / "figures"

rows = list(csv.DictReader(open(data_csv)))
years = [int(r["year"]) for r in rows]
trees = [int(r["trees_cumulative"]) for r in rows]
dated = [int(r["dated_trees_cumulative"]) for r in rows]
parts = [int(r["partitions_covered"]) for r in rows]
species = [int(r["species_in_canonical_trees"]) for r in rows]

mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False

fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)

# ── Left panel: cumulative trees + dated trees ──
ax = axes[0]
ax.plot(years, trees, "-o", color="#2a6fbf", markersize=4, linewidth=1.6, label="All trees")
ax.plot(years, dated, "-o", color="#1f6b29", markersize=4, linewidth=1.6, label="Time-calibrated")
ax.fill_between(years, 0, trees, alpha=0.08, color="#2a6fbf")
ax.set_xlabel("Year")
ax.set_ylabel("Cumulative atlas trees")
ax.set_xlim(min(years) - 0.5, max(years) + 0.5)
ax.set_xticks(range(min(years), max(years) + 1, 2))
ax.legend(loc="upper left", frameon=False, fontsize=9)
ax.set_title("a. Tree availability", loc="left", fontsize=11, fontweight="bold")

# ── Right panel: partitions covered + species in trees ──
ax = axes[1]
ax2 = ax.twinx()
ax2.spines["top"].set_visible(False)

line1 = ax.plot(years, parts, "-o", color="#8a4d0a", markersize=4, linewidth=1.6, label="Partitions with a tree")
line2 = ax2.plot(years, [s / 1000 for s in species], "-s", color="#1a4d8a",
                 markersize=4, linewidth=1.6, label="Species in canonical trees (thousands)")

ax.set_xlabel("Year")
ax.set_ylabel("Partitions covered")
ax2.set_ylabel("Species in canonical trees (×1000)")
ax.set_xlim(min(years) - 0.5, max(years) + 0.5)
ax.set_xticks(range(min(years), max(years) + 1, 2))
ax.set_title("b. Coverage of the tree of life", loc="left", fontsize=11, fontweight="bold")

# Combined legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc="upper left", frameon=False, fontsize=9)

# Save
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / "figure_atlas_growth.png", dpi=300, bbox_inches="tight")
plt.savefig(outdir / "figure_atlas_growth.pdf", bbox_inches="tight")
print(f"Wrote {outdir / 'figure_atlas_growth.png'}")
print(f"Wrote {outdir / 'figure_atlas_growth.pdf'}")
