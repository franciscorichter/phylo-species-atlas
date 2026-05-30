"""Figure S2: Data source landscape.

Two-panel figure showing the distribution of source trees by:
  (a) data repository
  (b) publication year

Data: Data/data_provenance.csv (the 49 source datasets, filtered to the
46 with deposit records).

Output: Paper/figures/figure_data_landscape.{png,pdf}
"""
import csv
import pathlib
import collections
import matplotlib.pyplot as plt
import matplotlib as mpl

script_dir = pathlib.Path(__file__).resolve().parent
root = script_dir.parent.parent
prov_csv = root.parent / "Data" / "data_provenance.csv"
outdir = root / "figures"

mpl.rcParams["font.family"] = "Helvetica"

with open(prov_csv) as f:
    rows = list(csv.DictReader(f))

# Filter to rows with deposit records (data_source non-empty) and a year
deposited = [
    r for r in rows
    if r.get("data_source", "").strip() and r.get("year", "").strip().isdigit()
]
n = len(deposited)

repos = collections.Counter(r["data_source"].strip() for r in deposited)
years = collections.Counter(int(r["year"]) for r in deposited)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.1, 1.4]})

# (a) Repository distribution — horizontal bars
ax = axes[0]
labels, counts = zip(*sorted(repos.items(), key=lambda x: -x[1]))
ypos = range(len(labels))
ax.barh(ypos, counts, color="#2C7FB8", edgecolor="white", height=0.7)
ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Number of datasets", fontsize=10)
ax.set_title(f"(a) Repository distribution ({n} datasets)", fontsize=11, loc="left", weight="bold")
for i, c in enumerate(counts):
    ax.text(c + 0.15, i, str(c), va="center", fontsize=9, color="#333")
ax.set_axisbelow(True)
ax.grid(axis="x", alpha=0.3, linewidth=0.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# (b) Publication-year distribution
ax = axes[1]
ymin, ymax = min(years), max(years)
yspan = list(range(ymin, ymax + 1))
counts_y = [years.get(y, 0) for y in yspan]
colors = ["#2C7FB8" if y >= 2018 else "#999999" for y in yspan]
ax.bar(yspan, counts_y, color=colors, edgecolor="white", width=0.7)
ax.set_xlabel("Publication year", fontsize=10)
ax.set_ylabel("Number of datasets", fontsize=10)
ax.set_title(f"(b) Publication year ({ymin}–{ymax})", fontsize=11, loc="left", weight="bold")
ax.set_xticks(yspan[::2])
ax.tick_params(axis="x", labelsize=9, rotation=0)
ax.tick_params(axis="y", labelsize=9)
ax.set_axisbelow(True)
ax.grid(axis="y", alpha=0.3, linewidth=0.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# annotation: "2018 onward"
post_2018 = sum(years[y] for y in years if y >= 2018)
ax.text(0.98, 0.95, f"{post_2018}/{n} datasets\nfrom 2018 onward",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9, color="#2C7FB8", weight="bold",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#2C7FB8", lw=0.8))

plt.tight_layout()
outdir.mkdir(parents=True, exist_ok=True)
plt.savefig(outdir / "figure_data_landscape.png", dpi=300, bbox_inches="tight")
plt.savefig(outdir / "figure_data_landscape.pdf", bbox_inches="tight")
print(f"Wrote {outdir / 'figure_data_landscape.png'}")
print(f"  ({n} deposited datasets across {len(repos)} repositories, {ymin}–{ymax})")
