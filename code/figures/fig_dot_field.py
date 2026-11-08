"""Figure 1 (dot-field version): described eukaryotic species at the four
levels of usability, drawn as the Species Phylogenies Catalogue draws them.

(a)-(d) The platform's sunflower disc of every described eukaryotic species,
        one dot per SCALE species (SCALE read from the platform's app.js),
        dots allotted to the groups by largest remainder and laid out in
        tree-of-life order (data/tree_of_groups.json).  Lit in brass: the
        species named in a reference tree (C0), placed with molecular data
        (C1), in a chronogram in absolute time (C1d), with date uncertainty
        archived (C1u).  Inside each group's sector the lit dots take the
        points nearest the centre, so a group's lit share is one region.
        Under each disc: percentage, count, and the loss from the previous
        level split by reason.
(e)     The platform's final view: nested discs at the density of (a)-(d).
        Named (brass) at the centre of the described (ink), faint dots out to
        the low estimate, a thin circle at the central and a dotted circle at
        the high sum of the groups' estimates.

Inputs (all read, none typed):
  Data/coverage_basis.csv                 described and level counts per group
  Data/data_estimates.csv                 estimated_low/total/high per group
  Paper/inputs/fig1_coverage_cascade.csv  losses between levels, by reason
  Paper/inputs/headline_numbers.json      cross-check of the totals
  release/files/site/app.js               SCALE (species per dot)
  release/files/site/data/tree_of_groups.json   order of the groups
Output: Paper/figures/figure_dot_field.{png,pdf}
"""
import csv
import json
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=6.5)

INK, BRASS, SOFT, FAINT, LINE = "#1f2328", "#b08d4f", "#8b9098", "#a3a8b0", "#e3e5e8"
SITE = next(p for p in (afc.REPO_ROOT / "release" / "files" / "site",
                        afc.REPO_ROOT / "site") if p.is_dir())

# ── the platform's constants and order ─────────────────────────────────────
m = re.search(r"const\s+SCALE\s*=\s*(\d+)", (SITE / "app.js").read_text())
SCALE = int(m.group(1))
tree = json.loads((SITE / "data" / "tree_of_groups.json").read_text())
order = []


def walk(n):
    if n.get("g"):
        order.append(n["g"])
    for c in n.get("c", []):
        walk(c)


walk(tree["tree"])

# ── the deposit ────────────────────────────────────────────────────────────
basis = afc.eukaryotic(afc.load_basis())
by = {r["partition"]: r for r in basis}
groups = [g for g in order if g in by]
if set(groups) != set(by):
    raise SystemExit("tree_of_groups.json and coverage_basis.csv disagree on "
                     "the eukaryotic groups")
denom = afc.denominator(basis)
tot = afc.rung_totals(basis)
head = afc.load_headlines()
afc.check_against_headlines(
    head, numerators=tot, denominator_euk=denom,
    percents={k: 100.0 * v / denom for k, v in tot.items()},
    label="fig_dot_field")
n_groups = len(basis)
n_contrib = {k: sum(1 for r in basis if r[k] > 0) for k in afc.RUNGS}

loss = {}
with open(afc.INPUTS_DIR / "fig1_coverage_cascade.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if r["kind"] == "loss":
            loss[r["group"]] = int(float(r["species"]))
STEP = {"C1": [("TaxPlaced", "named, not in the\nmolecular analysis")],
        "C1d": [("RelTime", "tree in relative time"),
                ("NotChron", "tree not a chronogram")],
        "C1u": [("NoInterval", "no date uncertainty archived"),
                ("Placement", "set varies tip placement only")]}
prev = {"C1": "C0", "C1d": "C1", "C1u": "C1d"}
for k, parts in STEP.items():
    if sum(loss[c] for c, _ in parts) != tot[prev[k]] - tot[k]:
        raise SystemExit(f"{k}: losses by reason do not add to the step loss")
if set(loss) != {c for v in STEP.values() for c, _ in v}:
    raise SystemExit("cascade losses do not match the reasons drawn")

est = {}
# the groups' estimates as the paper's Figure 3 data file carries them (the
# rebuild writes it from the curated per-group records)
with open(afc.INPUTS_DIR / "fig1_coverage_bars.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        est[r["group"]] = r


def est_sum(col):
    # as the platform: a group without an estimate counts its described species
    return sum(int(est[g][col]) if g in est and est[g][col] else by[g]["described"]
               for g in groups)


EST = {c: est_sum("estimated_" + c) for c in ("low", "total", "high")}


# ── dots per group, largest remainder, as app.js ───────────────────────────
def largest_remainder(values):
    target = round(sum(values) / SCALE)
    raw = [v / SCALE for v in values]
    base = [math.floor(x) for x in raw]
    rest = target - sum(base)
    for _, i in sorted(((raw[i] - base[i], i) for i in range(len(raw))),
                       key=lambda t: -t[0])[:rest]:
        base[i] += 1
    return base


dots = {g: {} for g in groups}
for g, n in zip(groups, largest_remainder([by[g]["described"] for g in groups])):
    dots[g]["all"] = n
p = "all"
for k in afc.RUNGS:
    for g, n in zip(groups, largest_remainder([by[g][k] for g in groups])):
        dots[g][k] = min(n, dots[g][p])
    p = k
n_desc = sum(d["all"] for d in dots.values())
n_named = sum(d["C0"] for d in dots.values())
low_dots = round(EST["low"] / SCALE)

# ── geometry (millimetres) ─────────────────────────────────────────────────
GOLD = math.pi * (3 - math.sqrt(5))
R = 13.5                                   # radius of the discs (a)-(d)
AREA = math.pi * R * R / n_desc            # mm^2 per dot, kept in (e)
rad = lambda nd: math.sqrt(nd * AREA / math.pi)
PITCH = math.sqrt(AREA)
OPEN = 0.05 * 2 * math.pi
START = -math.pi / 2 + OPEN / 2            # app.js: the ring starts at the top


def field():
    """Sunflower disc, points sorted by angle from START (app.js fieldTargets),
    then within each group's run of points the lit dots take the innermost."""
    i = np.arange(n_desc)
    r = R * np.sqrt((i + 0.5) / n_desc)
    a = i * GOLD
    x, y = np.cos(a) * r, np.sin(a) * r
    ang = np.mod(a - START, 2 * math.pi)
    o = np.argsort(ang, kind="stable")
    x, y, r = x[o], y[o], r[o]
    rank = {}                              # group -> indices, innermost first
    s = 0
    for g in groups:
        idx = np.arange(s, s + dots[g]["all"])
        rank[g] = idx[np.argsort(r[idx], kind="stable")]
        s += dots[g]["all"]
    return x, -y, rank                     # y flipped: matplotlib y is up


FX, FY, RANK = field()


PREV = {"C0": None, "C1": "C0", "C1d": "C1", "C1u": "C1d"}


def disc(ax, cx, cy, level):
    """Lit dots in brass; the dots lit at the previous level and lost at this
    one in pale brass; the rest of the described species in pale ink."""
    def mask(lv):
        m = np.zeros(n_desc, bool)
        if lv:
            for g in groups:
                m[RANK[g][:dots[g][lv]]] = True
        return m
    lit, was = mask(level), mask(PREV[level])
    lost = was & ~lit
    rest = ~lit & ~lost
    sz = (PITCH * 0.95 * 72 / 25.4) ** 2   # dot diameter 0.95 pitch, in pt^2
    ax.scatter(cx + FX[rest], cy + FY[rest], s=sz, c=INK, alpha=afc.UNLIT_ALPHA,
               linewidths=0)
    ax.scatter(cx + FX[lost], cy + FY[lost], s=sz, c=BRASS, alpha=0.6,
               linewidths=0)
    ax.scatter(cx + FX[lit], cy + FY[lit], s=sz * 1.1, c=BRASS,
               linewidths=0)


def pct(v):
    return f"{100 * v / denom:.1f}%"


W, H = 180.0, 127.5
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect("equal")
ax.axis("off")

LEVEL_NAME = {"C0": "Named in a reference tree", "C1": "Placed with molecular data",
              "C1d": "Dated in absolute time", "C1u": "Date uncertainty archived"}
cells = {"C0": (23.0, 97.0), "C1": (66.0, 97.0),
         "C1d": (23.0, 39.0), "C1u": (66.0, 39.0)}
for letter, k in zip("abcd", afc.RUNGS):
    cx, cy = cells[k]
    disc(ax, cx, cy, k)
    ax.text(cx - R - 1, cy + R + 4.5, f"({letter})", fontsize=7.5,
            fontweight="bold", color=INK, va="bottom")
    ax.text(cx - R + 5, cy + R + 4.5, LEVEL_NAME[k], fontsize=7, color=INK,
            va="bottom")
    ax.text(cx, cy - R - 2.2, pct(tot[k]), fontsize=9, color=BRASS,
            ha="center", va="top", fontweight="bold")
    ax.text(cx, cy - R - 7.2,
            f"{int(tot[k]):,} species · {n_contrib[k]} of {n_groups} groups",
            fontsize=6.5, color=INK, ha="center", va="top")
    if k in STEP:
        lines = [f"−{loss[c]:,} {why}" for c, why in STEP[k]]
        ax.text(cx, cy - R - 10.6, "\n".join(lines), fontsize=6, color=INK,
                ha="center", va="top", linespacing=1.25)

# the groups that carry at least 5 % of a panel's lit dots, named at that panel:
# a leader from the centre of the group's lit dots to a label beside the disc
LABELLED = {}


def name_lit(cx, cy, level):
    picks = []
    n_lit = sum(dots[g][level] for g in groups)
    for g in groups:
        k = dots[g][level]
        if k and k / n_lit >= 0.05:
            ii = RANK[g][:k]
            picks.append((g, cx + FX[ii].mean(), cy + FY[ii].mean()))
    # neighbours whose lit dots sit within 3 mm share one leader and one label
    clusters = []
    for g, px, py in picks:
        if clusters and math.hypot(px - clusters[-1][1], py - clusters[-1][2]) < 3.0:
            c = clusters[-1]
            n = len(c[0])
            clusters[-1] = (c[0] + [g], (c[1] * n + px) / (n + 1),
                            (c[2] * n + py) / (n + 1))
        else:
            clusters.append(([g], px, py))
    for side in (1, -1):
        mine = sorted([q for q in clusters if (q[1] - cx) * side >= 0],
                      key=lambda q: -q[2])
        last = None
        for names, px, py in mine:
            h = 0.0
            y = py if last is None else min(py, last - 2.60 - h / 2)
            lx = cx + side * (R + 1.5)
            ax.plot([px, lx - side * 0.5], [py, y], color=INK, lw=0.3,
                    alpha=0.6)
            big = max(names, key=lambda g: dots[g][level])
            more = f" +{len(names) - 1}" if len(names) > 1 else ""
            LABELLED.setdefault(level, []).append((big, names))
            ax.text(lx, y, big + more, fontsize=6, color=INK,
                    ha="left" if side > 0 else "right", va="center",
                    linespacing=1.15)
            last = y - h / 2


for k in afc.RUNGS:
    name_lit(*cells[k], k)

ax.text(2, H - 1, f"One dot = {SCALE:,} described eukaryotic species.   "
        f"(a)–(d): brass, species at that level; pale brass, lost since the level before.",
        fontsize=6, color=INK, va="top")

# ── (e) described against the estimated total ──────────────────────────────
ex, ey = 133.5, 62.0
Rlow, Rtot, Rhigh = rad(low_dots), rad(EST["total"] / SCALE), rad(EST["high"] / SCALE)
i = np.arange(low_dots)
r = np.sqrt((i + 0.5) * AREA / math.pi)
a = i * GOLD
x, y = ex + np.cos(a) * r, ey + np.sin(a) * r
sz = (PITCH * 0.85 * 72 / 25.4) ** 2
ax.scatter(x[n_desc:], y[n_desc:], s=sz, c=FAINT, alpha=0.3, linewidths=0)
ax.scatter(x[n_named:n_desc], y[n_named:n_desc], s=sz, c=INK, alpha=0.28,
           linewidths=0)
ax.scatter(x[:n_named], y[:n_named], s=sz * 1.15, c=BRASS, linewidths=0)
th = np.linspace(0, 2 * math.pi, 721)
ax.plot(ex + Rtot * np.cos(th), ey + Rtot * np.sin(th), color=INK, lw=0.5)
ax.plot(ex + Rhigh * np.cos(th), ey + Rhigh * np.sin(th), color=INK, lw=0.7,
        ls=(0, (0.6, 2.20)), dash_capstyle="round")
n_pub = sum(1 for g in groups if g in est and est[g]["estimated_high"]
            and int(est[g]["estimated_high"]) > by[g]["described"])
top = ey + Rhigh
ax.text(ex - Rhigh, top + 12.5, "(e)", fontsize=7.5, fontweight="bold",
        color=INK, va="bottom")
ax.text(ex - Rhigh + 6, top + 12.5,
        "Described species against the estimated total, at the same scale",
        fontsize=7, color=INK, va="bottom")
ax.text(ex - Rhigh + 6, top + 11.6,
        f"Each estimate sums the {n_groups} groups' values; "
        f"a group with no published estimate counts its described species.",
        fontsize=6, color=INK, va="top", linespacing=1.25)


def mio(v):
    return f"{v / 1e6:.1f} million"


def share(v):
    return (f"{100 * v / EST['low']:.1f}% of the low,\n"
            f"{100 * v / EST['high']:.1f}% of the high\nestimate")


# labels in the empty annulus between the low disc and the central circle
lab_x = ex + Rlow + 0.3


def label(y_text, r_from, ang, head, body, col):
    px, py = ex + math.cos(ang) * r_from, ey + math.sin(ang) * r_from
    ax.plot([px, lab_x - 0.8], [py, y_text], color=SOFT, lw=0.35)
    ax.text(lab_x, y_text + 0.9, head, fontsize=6.5, color=col,
            fontweight="bold", va="bottom")
    ax.text(lab_x, y_text + 0.4, body, fontsize=6, color=INK, va="top",
            linespacing=1.20)


Rnamed, Rdesc = rad(n_named), rad(n_desc)
label(ey + 14.0, Rnamed * 0.7, math.radians(40), "Named",
      f"{int(tot['C0']):,}", BRASS)
label(ey + 1.0, Rdesc * 0.92, math.radians(8), "Described",
      f"{denom:,}", INK)
label(ey - 12.0, Rlow * 0.97, math.radians(-30), "Low estimate",
      f"{mio(EST['low'])}", SOFT)
# the two circles, labelled inside the central circle: central at the top,
# high at the bottom with a leader across the thin circle to the dotted one
for rr, txt, sign in ((Rtot, f"Central estimate, {mio(EST['total'])}", 1),
                      (Rhigh, f"High estimate, {mio(EST['high'])}", -1)):
    qy = ey + sign * (Rtot - 6.0)
    ax.plot([ex, ex], [qy + sign * 1.4, ey + sign * rr], color=SOFT, lw=0.35)
    ax.text(ex, qy, txt, fontsize=6, color=INK, ha="center", va="center")

for lv, items in LABELLED.items():
    for big, names in items:
        if len(names) > 1:
            print(f"{lv}: {big} +{len(names) - 1} = {', '.join(names)}")
afc.savefig(plt, afc.FIGURES_DIR, "figure_dot_field")
