"""Figure: the platform's ring view -- described species around the tree of the groups.

A print rendering of the ring of the Species Phylogenies Catalogue
(release/files/site/app.js, steps "Tree" and "1 named").  The geometry is the
site's own:

  * one dot per SCALE described species (SCALE is read from app.js);
  * dots are allocated to groups by the largest-remainder method, eukaryotes
    and prokaryotes separately, and a group's named dots never exceed its
    described dots (app.js loadData);
  * groups sit in the leaf order of the Open Tree of Life backbone in
    release/files/site/data/tree_of_groups.json, one arc per group, K dots
    deep, arc width proportional to described species, a gap at the top;
  * dots fill each arc row by row from the inside, so the named share is a
    band on the inner edge (app.js layoutRing);
  * the backbone is a radial cladogram inside the ring, tips at the middle of
    each arc, node radius proportional to depth (app.js layoutRing).

Encoding (one quantity in colour: named species):
  brass filled  -- described species named in the group's reference tree (C0)
  ink filled    -- described species of a group that has a reference tree, not named
  open grey     -- described species of a group with no species-level reference tree
  faint grey    -- Bacteria and Archaea, GTDB species clusters, outside every
                   eukaryotic percentage

Every number is read from Data/coverage_basis.csv and checked against
Paper/inputs/headline_numbers.json.

Output: Paper/figures/figure_tree_of_groups.{png,pdf}
"""
import json
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Circle

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

SITE = next(p for p in (afc.REPO_ROOT / "release" / "files" / "site",
                        afc.REPO_ROOT / "site") if p.is_dir())
APP_JS = SITE / "app.js"
TREE_JSON = SITE / "data" / "tree_of_groups.json"

# Platform palette (style.css / app.js)
INK = "#1f2328"
BRASS = "#b08d4f"
SOFT = "#8b9098"
FAINT = "#a3a8b0"
LINE = "#e3e5e8"

# ── Data ───────────────────────────────────────────────────────────────────
m = re.search(r"const\s+SCALE\s*=\s*(\d+)\s*;", APP_JS.read_text())
if not m:
    raise SystemExit("SCALE not found in app.js")
SCALE = int(m.group(1))

basis = afc.load_basis()
head = afc.load_headlines()
euk = afc.eukaryotic(basis)
prok = afc.prokaryotic(basis)
denom = afc.denominator(euk)
afc.check_against_headlines(head, numerators={"C0": afc.rung_totals(euk)["C0"]},
                            percents={"C0": afc.coverage_pct(euk)["C0"]},
                            denominator_euk=denom, label="fig_tree_of_groups")

tog = json.loads(TREE_JSON.read_text())
TREE, BINS = tog["tree"], tog["bins"]
G = {r["partition"]: r for r in basis}


def largest_remainder(values, scale):
    target = round(sum(values) / scale)
    raw = [v / scale for v in values]
    base = [math.floor(x) for x in raw]
    rest = target - sum(base)
    order_ = sorted(range(len(raw)), key=lambda i: -(raw[i] - base[i]))
    for i in order_[:rest]:
        base[i] += 1
    return base


dots = {}
for part in (euk, prok):
    alln = largest_remainder([r["described"] for r in part], SCALE)
    named = largest_remainder([r["C0"] for r in part], SCALE)
    for r, a, n in zip(part, alln, named):
        dots[r["partition"]] = {"all": a, "named": min(n, a)}

order = []


def walk(n):
    if "g" in n:
        order.append(n["g"])
    for c in n.get("c", []):
        walk(c)


walk(TREE)
missing = set(order) ^ set(G)
if missing:
    raise SystemExit(f"tree_of_groups.json and coverage_basis.csv disagree: {missing}")

# ── Geometry (app.js layoutRing, in mm) ────────────────────────────────────
FIG_W, FIG_H = 180.0, 172.0
cx, cy = FIG_W / 2 + 8, FIG_H / 2 + 1
Ro = 55.0
OPEN = 0.05 * 2 * math.pi
GAP_COLS = 0.7
best = None
for K in range(4, 13):
    cols = sum(max(1, math.ceil(dots[n]["all"] / K)) for n in order) + GAP_COLS * len(order)
    p = 0.3 * Ro / K
    for _ in range(30):
        Rm = Ro - K * p / 2
        p = min(0.3 * Ro / K, (2 * math.pi - OPEN) * Rm / cols)
    if best is None or p > best[1]:
        best = (K, p, cols)
K, p, cols = best
Ri = Ro - K * p
dth = (2 * math.pi - OPEN) / cols
start = -math.pi / 2 + OPEN / 2
dotR = max(0.2, p * 0.36)


def xy(a, r):
    # canvas angles run clockwise with y down; mirror y for matplotlib
    return cx + math.cos(a) * r, cy - math.sin(a) * r


arcs, col = {}, 0.0
cats = {"named": [], "ink": [], "open": [], "prok": []}
for n in order:
    d, row = dots[n], G[n]
    c = max(1, math.ceil(d["all"] / K))
    a0 = start + col * dth
    arcs[n] = {"a0": a0, "a1": a0 + c * dth, "mid": a0 + c * dth / 2, "cols": c}
    for k in range(d["all"]):
        a = a0 + (k % c + 0.5) * dth
        r = Ri + (k // c + 0.5) * p
        if row["is_prokaryotic"]:
            cats["prok"].append(xy(a, r))
        elif not row["has_tree"]:
            cats["open"].append(xy(a, r))
        elif k < d["named"]:
            cats["named"].append(xy(a, r))
        else:
            cats["ink"].append(xy(a, r))
    col += c + GAP_COLS

# backbone: radial cladogram
Rt = Ri - max(2.0, p * 1.5)


def height(n):
    return 0 if "g" in n else 1 + max(height(c) for c in n["c"])


H_TREE = height(TREE)
segs = []


def place(n, depth):
    if "g" in n:
        n["_ang"], n["_rad"] = arcs[n["g"]]["mid"], Rt
        return
    for c in n["c"]:
        place(c, depth + 1)
    angs = [c["_ang"] for c in n["c"]]
    n["_ang"] = (min(angs) + max(angs)) / 2
    n["_rad"] = Rt * depth / H_TREE * 0.97


place(TREE, 0)


def draw_tree(n):
    if "g" in n:
        return
    a0 = min(c["_ang"] for c in n["c"])
    a1 = max(c["_ang"] for c in n["c"])
    steps = max(2, int(abs(a1 - a0) / 0.02))
    segs.append([xy(a0 + (a1 - a0) * i / steps, n["_rad"]) for i in range(steps + 1)])
    for c in n["c"]:
        segs.append([xy(c["_ang"], n["_rad"]), xy(c["_ang"], c["_rad"])])
        draw_tree(c)


draw_tree(TREE)

# ── Draw ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W / 25.4, FIG_H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")

ax.add_collection(LineCollection(segs, colors=INK, linewidths=0.35, alpha=0.4,
                                 capstyle="butt"))


def dotset(pts, **kw):
    ax.add_collection(PatchCollection([Circle(q, dotR) for q in pts], **kw))


dotset(cats["ink"], facecolor=INK, edgecolor="none", alpha=afc.UNLIT_ALPHA)
dotset(cats["named"], facecolor=BRASS, edgecolor="none")
dotset(cats["open"], facecolor="white", edgecolor=INK, linewidth=0.35, alpha=0.8)
dotset(cats["prok"], facecolor=LINE, edgecolor="none")

# groups smaller than one dot: a ring at the tip (app.js drawBackbone)
for n in order:
    if dots[n]["all"] == 0:
        q = xy(arcs[n]["mid"], Ri + p * 0.6)
        ec = INK if not G[n]["has_tree"] else INK
        ax.add_patch(Circle(q, dotR * 1.25, facecolor="white", edgecolor=ec, linewidth=0.4))

# labels: radial, just outside each arc (app.js drawLabels places them at the
# arc's middle).  Every group with a reference tree is labelled, and a group
# without one only if it holds at least MIN_LABEL_DOTS dots; the rest are
# listed by the script for the caption.  Labels are relaxed along the ring so
# neighbours never overlap; a leader line is drawn only for a label pushed
# away from its arc.
LABEL_PT = 6.3
MIN_LABEL_DOTS = 10          # = 5,000 described species at the site's scale
LAB_NEAR, LAB_FAR = Ro + 1.20, Ro + 4.2
MIN_SEP_MM = 2.60
SHIFT_MM = 0.9               # displacement beyond which a leader is drawn


def label_text(n):
    row = G[n]
    if row["is_prokaryotic"] or not row["has_tree"]:
        t, color = n, SOFT
    else:
        pct = 100.0 * row["C0"] / row["described"]
        t, color = (f"{n}  {pct:.0f}%" if pct >= 1 else f"{n}  <1%"), INK
    if n in BINS:
        t += "*"
    return t, color


labelled = [n for n in order if G[n]["has_tree"] or dots[n]["all"] >= MIN_LABEL_DOTS]
unlabelled = [n for n in order if n not in labelled]
want = [arcs[n]["mid"] for n in labelled]
pos = want[:]
sep = MIN_SEP_MM / LAB_NEAR
for _ in range(3000):
    moved = False
    for i in range(1, len(pos)):
        gap = pos[i] - pos[i - 1]
        if gap < sep:
            sh = (sep - gap) / 2 + 1e-6
            pos[i - 1] -= sh
            pos[i] += sh
            moved = True
    pos = [q + 0.02 * (a - q) for q, a in zip(pos, want)]
    if not moved:
        break

for n, a_arc, a_lab in zip(labelled, want, pos):
    shifted = abs(a_lab - a_arc) * LAB_NEAR > SHIFT_MM
    r_lab = LAB_FAR if shifted else LAB_NEAR
    if shifted:
        p0, p1 = xy(a_arc, Ro + 0.6), xy(a_arc, Ro + 1.8)
        p2 = xy(a_lab, r_lab - 0.5)
        ax.plot([p0[0], p1[0], p2[0]], [p0[1], p1[1], p2[1]], color=FAINT,
                lw=0.3, solid_capstyle="butt")
    x, y = xy(a_lab, r_lab)
    deg = -math.degrees(a_lab)
    right = math.cos(a_lab) >= 0
    t, color = label_text(n)
    ax.text(x, y, t, fontsize=LABEL_PT, color=color,
            rotation=deg if right else deg + 180, rotation_mode="anchor",
            ha="left" if right else "right", va="center")

# ── Clade arcs: a thin arc just inside the ring, spanning a clade's groups ──
# Membership is read, never chosen here:
#   * a coverage_basis category, as the maximal run of consecutive groups of
#     that category in ring order (Plants, Vertebrates, Arthropods, and the
#     protist run of "Microbes & protists" between the plants and the animals);
#   * a node of tree_of_groups.json, as the leaves of the smallest node that
#     holds every member of a catch-all bin ("Other Mollusca" -> molluscs).
def category_runs():
    runs, cur = [], []
    for n in order:
        if cur and G[n]["category"] == G[cur[-1]]["category"]:
            cur.append(n)
        else:
            if cur:
                runs.append(cur)
            cur = [n]
    runs.append(cur)
    return runs


def leaves(n):
    return [n["g"]] if "g" in n else [x for c in n["c"] for x in leaves(c)]


def smallest_node_with(members, n=TREE):
    for c in n.get("c", []):
        if set(members) <= set(leaves(c)):
            return smallest_node_with(members, c)
    return n


CLADE_FROM_CATEGORY = {"Plants": "Plants", "Vertebrates": "Vertebrates",
                       "Arthropods": "Arthropods", "Microbes & protists": "Protists"}
clades = []
for run in category_runs():
    cat = G[run[0]]["category"]
    if cat in CLADE_FROM_CATEGORY and len(run) >= 3 and not any(G[n]["is_prokaryotic"] for n in run):
        clades.append((CLADE_FROM_CATEGORY[cat], run, "category"))
moll = leaves(smallest_node_with(BINS["Other Mollusca"]))
clades.append(("Molluscs", [n for n in order if n in moll], "tree_of_groups.json"))

R_CLADE = Ri - 1.0
for name, run, src in clades:
    a0, a1 = arcs[run[0]]["a0"] + 0.15 * dth, arcs[run[-1]]["a1"] - 0.15 * dth
    k = max(2, int((a1 - a0) / 0.01))
    pts = [xy(a0 + (a1 - a0) * i / k, R_CLADE) for i in range(k + 1)]
    ax.plot([q[0] for q in pts], [q[1] for q in pts], color=SOFT, lw=0.6,
            solid_capstyle="butt", zorder=3)
    am = (a0 + a1) / 2
    x, y = xy(am, R_CLADE - 2.0)
    deg = -math.degrees(am) - 90          # tangential
    if math.sin(am) > 0:                  # lower half (canvas y down): keep upright
        deg += 180
    ax.text(x, y, name, fontsize=6.0, color=SOFT, rotation=deg,
            rotation_mode="anchor", ha="center", va="center", zorder=4,
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none"))
    print(f"clade {name} ({src}): {', '.join(run)}")

# ── Numbers, all computed ──────────────────────────────────────────────────
c0 = sum(r["C0"] for r in euk)
n_euk = len(euk)
n_tree = sum(1 for r in euk if r["has_tree"])
n_none = n_euk - n_tree
desc_none = sum(r["described"] for r in euk if not r["has_tree"])
core = [r for r in euk if r["category"] == "Vertebrates" or r["partition"] == "Seed plants"]
rest = [r for r in euk if r not in core]
core_pct = 100.0 * sum(r["C0"] for r in core) / afc.denominator(core)
rest_pct = 100.0 * sum(r["C0"] for r in rest) / afc.denominator(rest)
prok_desc = afc.denominator(prok)
n_euk_dots = sum(dots[r["partition"]]["all"] for r in euk)

TXT = 7.0
core_share = 100.0 * afc.denominator(core) / denom
ax.text(4, FIG_H - 4,
        f"Named in their group's reference tree: {100 * c0 / denom:.1f}% of {denom:,} described eukaryotic species\n"
        f"{core_pct:.0f}% of the species of vertebrates and seed plants ({core_share:.0f}% of described species)\n"
        f"{rest_pct:.1f}% of the species of all other eukaryotes ({100 - core_share:.0f}% of described species)\n"
        f"{n_none} of {n_euk} groups have no species-level tree (open dots, grey labels)",
        fontsize=TXT + 0.5, color=INK, va="top", ha="left", linespacing=1.4)

# key, bottom left
dy = 3.5
key = [
    (dict(facecolor=BRASS, edgecolor="none"), "named in the group's reference tree"),
    (dict(facecolor=INK, edgecolor="none", alpha=afc.UNLIT_ALPHA), "not named; the group has a reference tree"),
    (dict(facecolor="white", edgecolor=INK, linewidth=0.45, alpha=0.8),
     f"the group has no species-level tree ({n_none} of {n_euk} groups)"),
    (dict(facecolor=LINE, edgecolor="none"),
     f"Bacteria, Archaea: {prok_desc:,} GTDB clusters, outside every %"),
]
lines = [f"one dot = {SCALE:,} described species; one arc = one group,",
         "placed by the tree of the groups drawn inside the ring",
         "* catch-all group, not a clade: at the common ancestor of its members"]
kx = 5.6
y0 = 4.0 + (len(key) + len(lines) - 1) * dy
for i, t in enumerate(lines[:2]):
    ax.text(4, y0 - i * dy, t, fontsize=TXT, color=INK, va="center", ha="left")
for i, (kw, lab) in enumerate(key):
    yy = y0 - (i + 2) * dy
    ax.add_patch(Circle((kx, yy), 0.9, **kw))
    ax.text(kx + 2.4, yy, lab, fontsize=TXT, color=INK, va="center", ha="left")
ax.text(4, y0 - (len(key) + 2) * dy, lines[2], fontsize=TXT, color=INK, va="center", ha="left")

afc.savefig(plt, afc.FIGURES_DIR, "figure_tree_of_groups", facecolor="white")

print(f"SCALE={SCALE} K={K} pitch={p:.3f}mm eukaryotic dots={n_euk_dots}")
print(f"C0={c0:,.0f} denom={denom:,} pct={100 * c0 / denom:.2f}")
print(f"groups: {n_euk} eukaryotic, {n_tree} with tree, {n_none} without; "
      f"described without tree {desc_none:,} ({100 * desc_none / denom:.1f}%)")
print(f"vert+seed plants named {core_pct:.2f}%  rest {rest_pct:.2f}%  prok {prok_desc:,}")
print("bins:", ", ".join(BINS))
print(f"labelled {len(labelled)}; unlabelled (MIN_LABEL_DOTS={MIN_LABEL_DOTS}), clockwise from the top:", ", ".join(unlabelled))
