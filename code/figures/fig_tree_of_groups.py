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
    each arc (app.js layoutRing).

What the print rendering changes, so that the picture is read without its
caption:
  * the cladogram steps inward by a node's height above its tips (one step
    per level), with heavier lines for the deeper splits, so a clade is the
    fan below its node; the site steps outward from the root;
  * five branches are named in bold inside the ring, each with the share of
    its described species that are named, a dot on its node and a thin arc
    along its groups.  Membership is read from the category column of
    coverage_basis.csv and from tree_of_groups.json, never typed;
  * every label is horizontal, one size, and is placed only where its measured
    extent clears the ring, the tree and every other label.  A group with a
    tree whose own arc has no room is labelled beside it with a short leader.
    A run of groups with trees too narrow for that (the vertebrates) gets one
    bracket on its arc and an aligned list beside the ring, in ring order; a
    run of unlabelled groups with no tree gets a grey bracket and their
    number; Bacteria and Archaea get one bracket and a label in ink;
  * a group smaller than one dot is a small mark at its tip: filled if it has
    a reference tree, open if it has none (the site draws a ring for both);
  * the dots are slightly larger for the same pitch (DOT_FILL).
The script prints every text it draws and the groups left without a label.

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
import numpy as np

import atlas_figure_common as afc

mpl.use("Agg")
afc.apply_style(mpl, size=7)

SITE = next(p for p in (afc.REPO_ROOT / "release" / "files" / "site",
                        afc.REPO_ROOT / "site") if p.is_dir())
APP_JS = SITE / "app.js"
TREE_JSON = SITE / "data" / "tree_of_groups.json"

# Platform palette (style.css / app.js), held in atlas_figure_common
INK, BRASS, SOFT, LINE = afc.INK, afc.BRASS, afc.SOFT, afc.LINE

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
FIG_W, FIG_H = 177.0, 146.0
Ro = 62.0
cx, cy = 87.25, Ro + 10.5
OPEN = 0.05 * 2 * math.pi
GAP_COLS = 0.7
DOT_FILL = 0.40            # dot radius over pitch (the site: 0.36)
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
dotR = max(0.2, p * DOT_FILL)


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

# ── Backbone: radial cladogram, tips on the ring, one step inward per level ─
Rt = Ri - max(1.75, p * 1.5)
R_ROOT = 3.25


def place(n):
    if "g" in n:
        n["_h"], n["_ang"], n["_rad"] = 0, arcs[n["g"]]["mid"], Rt
        return
    for c in n["c"]:
        place(c)
    angs = [c["_ang"] for c in n["c"]]
    n["_h"] = 1 + max(c["_h"] for c in n["c"])
    n["_ang"] = (min(angs) + max(angs)) / 2


place(TREE)
H_TREE = TREE["_h"]
STEP = (Rt - R_ROOT) / H_TREE


def set_rad(n):
    if "g" not in n:
        n["_rad"] = Rt - STEP * n["_h"]
        for c in n["c"]:
            set_rad(c)


set_rad(TREE)

TREE_LW_TIP, TREE_LW_ROOT = 0.28, 0.80


def tree_lw(h):
    return TREE_LW_TIP + (TREE_LW_ROOT - TREE_LW_TIP) * h / H_TREE


segs, seg_lw = [], []


def draw_tree(n):
    if "g" in n:
        return
    a0 = min(c["_ang"] for c in n["c"])
    a1 = max(c["_ang"] for c in n["c"])
    steps = max(2, int(abs(a1 - a0) / 0.02))
    segs.append([xy(a0 + (a1 - a0) * i / steps, n["_rad"]) for i in range(steps + 1)])
    seg_lw.append(tree_lw(n["_h"]))
    for c in n["c"]:
        segs.append([xy(c["_ang"], n["_rad"]), xy(c["_ang"], c["_rad"])])
        seg_lw.append(tree_lw(c["_h"]))
        draw_tree(c)


draw_tree(TREE)

# points along every tree line, for the label-clearance test
tree_pts = []
for line in segs:
    for (xa, ya), (xb, yb) in zip(line[:-1], line[1:]):
        k = max(1, int(math.hypot(xb - xa, yb - ya) / 0.3))
        tree_pts.extend((xa + (xb - xa) * i / k, ya + (yb - ya) * i / k) for i in range(k + 1))

# ── Clades named in the tree ───────────────────────────────────────────────
# Membership is read, never chosen here:
#   * a coverage_basis category, as the maximal run of consecutive groups of
#     that category in ring order (Plants, Vertebrates, Arthropods, and the
#     protist run of "Microbes & protists" between the plants and the animals).
#     The node of tree_of_groups.json whose leaves are that run is marked;
#     where the run is not one clade of the backbone (the protists) no node
#     is marked and the name is set by the largest clade inside the run;
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


def nodes(n):
    yield n
    for c in n.get("c", []):
        yield from nodes(c)


def smallest_node_with(members, n=TREE):
    for c in n.get("c", []):
        if set(members) <= set(leaves(c)):
            return smallest_node_with(members, c)
    return n


def largest_node_within(members):
    inside = [n for n in nodes(TREE) if "g" not in n and set(leaves(n)) <= set(members)]
    return max(inside, key=lambda n: len(leaves(n)))


def share(names):
    rows = [G[n] for n in names]
    return 100.0 * sum(r["C0"] for r in rows) / sum(r["described"] for r in rows)


def pct_text(v):
    return f"{v:.0f}%" if v >= 1 else "<1%"


RUNS = category_runs()
CLADE_FROM_CATEGORY = {"Plants": "Plants", "Vertebrates": "Vertebrates",
                       "Arthropods": "Arthropods", "Microbes & protists": "Protists"}
clades = []
for run in RUNS:
    cat = G[run[0]]["category"]
    if cat in CLADE_FROM_CATEGORY and len(run) >= 3 and not any(G[n]["is_prokaryotic"] for n in run):
        node = smallest_node_with(run)
        is_clade = set(leaves(node)) == set(run)
        src = "category run, one clade of the backbone"
        if not is_clade:
            node = largest_node_within(run)
            src = "category run, not one clade: no node marked, name set by its largest clade"
        if set(run) == {r["partition"] for r in euk if r["category"] == cat}:
            want = head["per_category_eukaryotic"][cat]["C0_pct"]
            if not math.isclose(share(run), want, abs_tol=1e-6):
                raise SystemExit(f"{cat}: share {share(run)!r} != headline {want!r}")
        clades.append((CLADE_FROM_CATEGORY[cat], node, run, src, is_clade))
moll = smallest_node_with(BINS["Other Mollusca"])
clades.append(("Molluscs", moll, [n for n in order if n in leaves(moll)], "tree_of_groups.json", True))

# ── Canvas ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W / 25.4, FIG_H / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("white")
renderer = fig.canvas.get_renderer()

LABEL_PT = 7.5        # every label, list row, key line and (in bold) branch name
HEAD_PT = 8.0         # the one headline, reduced if needed to span the canvas
MIN_PT = 6.0          # no text below this at print size
texts = []            # (rect, string) of every text drawn, for the collision check


def offsets(s, **kw):
    """Extent of a horizontal text relative to its anchor, in mm."""
    t = ax.text(0, 0, s, **kw)
    bb = t.get_window_extent(renderer)
    (x0, y0), (x1, y1) = ax.transData.inverted().transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    t.remove()
    return x0, y0, x1, y1


def rect_at(x, y, off):
    return x + off[0], y + off[1], x + off[2], y + off[3]


def overlap(a, b, gap=0.0):
    return not (a[2] + gap <= b[0] or b[2] + gap <= a[0] or a[3] + gap <= b[1] or b[3] + gap <= a[1])


def put(x, y, s, **kw):
    off = offsets(s, **kw)
    ax.text(x, y, s, **kw)
    texts.append((rect_at(x, y, off), s))
    return texts[-1][0]


def ring_clear(rc, radius):
    """True if the rectangle lies wholly outside the circle of that radius."""
    qx = min(max(cx, rc[0]), rc[2])
    qy = min(max(cy, rc[1]), rc[3])
    return math.hypot(qx - cx, qy - cy) >= radius


# ── Marks ──────────────────────────────────────────────────────────────────
ax.add_collection(LineCollection(segs, colors=SOFT, linewidths=seg_lw, capstyle="round",
                                 joinstyle="round", zorder=1))


def dotset(pts, **kw):
    ax.add_collection(PatchCollection([Circle(q, dotR) for q in pts], zorder=2, **kw))


MARK = {
    "named": dict(facecolor=BRASS, edgecolor="none"),
    "ink": dict(facecolor=INK, edgecolor="none", alpha=afc.UNLIT_ALPHA),
    "open": dict(facecolor="white", edgecolor=INK, linewidth=0.35, alpha=0.85),
    "prok": dict(facecolor=LINE, edgecolor="none"),
}
for kind in ("ink", "named", "open", "prok"):
    dotset(cats[kind], **MARK[kind])

# a group smaller than one dot: a small mark at its tip, open if it has no
# tree, filled otherwise (brass if most of its species are named)
for n in order:
    if dots[n]["all"] == 0:
        row = G[n]
        q = xy(arcs[n]["mid"], Ri + p * 0.6)
        if not row["has_tree"]:
            kw = MARK["open"]
        elif 2 * row["C0"] >= row["described"]:
            kw = MARK["named"]
        else:
            kw = MARK["ink"]
        ax.add_patch(Circle(q, dotR * 0.8, zorder=2, **kw))

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

# ── Headline (one line) and key (the marks, the scale, the tree) ───────────
MARGIN = 1.5
headline = (f"{100 * c0 / denom:.1f}% of {denom:,} described eukaryotic species are named in "
            f"a species-level tree of their group; {n_none} of the {n_euk} groups have no such tree")
head_w = offsets(headline, fontsize=HEAD_PT, ha="left")[2]
head_pt = min(HEAD_PT, HEAD_PT * (FIG_W - 2 * MARGIN) / head_w)         # the one line spans the canvas at most
if min(head_pt, LABEL_PT) < MIN_PT:
    raise SystemExit(f"text below {MIN_PT} pt: headline {head_pt:.2f}, labels {LABEL_PT}")
reserved = [put(MARGIN, FIG_H - MARGIN, headline, fontsize=head_pt, color=INK,
                va="top", ha="left")]

KEY_Y1, KEY_Y0 = 6.75, 2.9
KEY_DOT, KEY_GAP = 0.85, 4.75
KW = dict(fontsize=LABEL_PT, color=INK, va="center", ha="left")
x = MARGIN
for kind, lab in (("named", "named in the group's species-level tree"),
                  ("ink", "not named in it"),
                  ("open", "the group has no such tree (name in grey)")):
    ax.add_patch(Circle((x + KEY_DOT, KEY_Y1), KEY_DOT, **MARK[kind]))
    x = put(x + 2 * KEY_DOT + 1.0, KEY_Y1, lab, **KW)[2] + KEY_GAP
x = put(MARGIN, KEY_Y0, f"One dot = {SCALE:,} described species.", **KW)[2] + KEY_GAP
ax.plot([x, x + 4.0], [KEY_Y0, KEY_Y0], color=SOFT, lw=tree_lw(H_TREE), solid_capstyle="butt")
x = put(x + 5.0, KEY_Y0, "tree of the groups: branching order only", **KW)[2] + KEY_GAP
ax.add_patch(Circle((x + 0.6, KEY_Y0), 0.6, facecolor=INK, edgecolor="none"))
put(x + 2.25, KEY_Y0, "a branch and, in bold, the share of its species named", **KW)
reserved += [r for r, _ in texts[1:]]

# ── Group labels ───────────────────────────────────────────────────────────
# Horizontal, one size.  A label sits just outside its own arc and is drawn
# only if its measured extent clears the ring and every label placed before
# it; larger groups are placed first.  A group with a tree that finds no room
# on its own arc is set beside it, at the end of a short leader.  A category
# run in which two or more groups with a tree cannot be labelled on their own
# arcs gets one bracket and an aligned list beside the ring, in ring order.
# A group with no tree is named (in grey) only if it holds MIN_LABEL_DOTS dots.
MIN_LABEL_DOTS = 25          # = 12,500 described species at the site's scale
LAB_R = Ro + 1.25            # clearance of a label from the outer edge of the ring
LAB_GAP = 0.95               # clearance between two labels
BRACKET_R = Ro + 1.15
LEADER_COLS = 3              # a label beside its arc sits at most this many columns away
LEADER_OUT = 1.25            # and that much further out, at the end of a short leader
ROW = 1.3 * (offsets("Ag", fontsize=LABEL_PT)[3] - offsets("Ag", fontsize=LABEL_PT)[1])


def label_text(n):
    row = G[n]
    if not row["has_tree"]:
        return n, SOFT
    return f"{n}  {pct_text(100.0 * row['C0'] / row['described'])}", INK


def bracket(a0, a1, color):
    """A thin arc just outside the ring over the arcs from a0 to a1, ends turned in."""
    k = max(2, int((a1 - a0) / 0.01))
    pts = [xy(a0, BRACKET_R - 0.7)] + [xy(a0 + (a1 - a0) * i / k, BRACKET_R) for i in range(k + 1)] \
        + [xy(a1, BRACKET_R - 0.7)]
    ax.plot([q[0] for q in pts], [q[1] for q in pts], color=color, lw=0.5,
            solid_capstyle="butt", solid_joinstyle="miter", zorder=3)


def direct_label(n, placed, beside=False):
    """Anchor and alignment of a horizontal label on the group's own arc, or
    (beside=True) within LEADER_COLS columns of it; None if neither is clear."""
    t, color = label_text(n)
    return arc_label(t, color, arcs[n]["a0"], arcs[n]["a1"], placed, beside)


def arc_label(t, color, a0, a1, placed, beside=False, r0=None):
    """Anchor, alignment and extent of a horizontal text just outside the arcs
    from a0 to a1, at least r0 from the centre; None if no anchor is clear."""
    mid, half = (a0 + a1) / 2, (a1 - a0) / 2
    r0 = LAB_R if r0 is None else r0
    if beside:
        shifts = [sg * (half + k * dth / 2) for k in range(1, 2 * LEADER_COLS + 1) for sg in (1, -1)]
    else:
        shifts = [f * half for f in (0, -0.25, 0.25, -0.5, 0.5, -0.75, 0.75, -1, 1)]
    for sh in shifts:
        a = mid + sh
        c, s = math.cos(a), math.sin(a)
        kw = dict(fontsize=LABEL_PT, color=color,
                  ha="left" if c > 0.4 else "right" if c < -0.4 else "center",
                  va="center" if abs(c) > 0.4 else "bottom" if s < 0 else "top")
        off = offsets(t, **kw)
        r = r0 + (LEADER_OUT if beside else 0)
        while not ring_clear(rect_at(*xy(a, r), off), r0 + (LEADER_OUT if beside else 0)):
            r += 0.1
        rc = rect_at(*xy(a, r), off)
        if rc[0] < MARGIN or rc[2] > FIG_W - MARGIN or rc[3] > FIG_H - MARGIN:
            continue
        if any(overlap(rc, q, LAB_GAP) for q in placed):
            continue
        return xy(a, r), t, kw, rc, (a if beside else None)
    return None


def priority(n):
    return (not G[n]["has_tree"], -G[n]["described"])


candidates = sorted((n for n in order if not G[n]["is_prokaryotic"]
                     and (G[n]["has_tree"] or dots[n]["all"] >= MIN_LABEL_DOTS)), key=priority)


def place_direct(skip, placed):
    out = {}
    todo = [n for n in candidates if n not in skip]
    with_tree = [n for n in todo if G[n]["has_tree"]]
    # groups with a tree first, on their arc or else beside it; then the others
    no_tree = [n for n in todo if not G[n]["has_tree"]]
    for group, beside in ((with_tree, False), (with_tree, True),
                          (no_tree, False), (no_tree, True)):
        for n in group:
            if n not in out:
                hit = direct_label(n, placed, beside)
                if hit:
                    out[n] = hit
                    placed.append(hit[3])
    return out


# Bacteria and Archaea: one bracket on their arcs and one label that says they
# are outside the count, in ink so that a grey name means only "no tree"
prok_names = [n for n in order if G[n]["is_prokaryotic"]]
a0, a1 = arcs[prok_names[0]]["a0"], arcs[prok_names[-1]]["a1"]
if any(not G[n]["is_prokaryotic"] for n in order[order.index(prok_names[0]):order.index(prok_names[-1]) + 1]):
    raise SystemExit("the prokaryotic groups are not one run of the ring")
PROK_LABEL = " and ".join(prok_names) + ": not counted"
hit = arc_label(PROK_LABEL, INK, a0, a1, list(reserved), r0=BRACKET_R + 0.9)
if hit is None:
    raise SystemExit("no clear place for the label of Bacteria and Archaea")
(x, y), t, kw, rc, _ = hit
ax.text(x, y, t, **kw)
texts.append((rc, t))
reserved.append(rc)
bracket(a0, a1, INK)

trial = place_direct(set(), list(reserved))
listed_runs = [run for run in RUNS if len(run) >= 3
               and sum(1 for n in run if G[n]["has_tree"] and not G[n]["is_prokaryotic"]
                       and n not in trial) >= 2]

in_list = set()
for run in listed_runs:
    a0, a1 = arcs[run[0]]["a0"], arcs[run[-1]]["a1"]
    am = (a0 + a1) / 2
    right = math.cos(am) >= 0
    rows = [n for n in run if G[n]["has_tree"]]
    if not right:
        rows = rows[::-1]                      # clockwise runs upward on the left
    in_list.update(run)
    pcts = [pct_text(100.0 * G[n]["C0"] / G[n]["described"]) for n in rows]
    kw = dict(fontsize=LABEL_PT, color=INK, va="center")
    wn = max(offsets(n, ha="left", **kw)[2] for n in rows)
    wp = max(-offsets(q, ha="right", **kw)[0] for q in pcts)
    sgn = 1 if right else -1
    x_in = cx + sgn * (BRACKET_R + 1.75)       # edge of the list nearest the ring
    x_out = x_in + sgn * (wn + 1.5 + wp)
    x_name, x_pct = (x_in, x_out) if right else (x_out, x_in)
    y_top = xy(am, Ro)[1] + ROW * (len(rows) - 1) / 2
    for i, (n, q) in enumerate(zip(rows, pcts)):
        y = y_top - i * ROW
        reserved.append(put(x_name, y, n, ha="left", **kw))
        reserved.append(put(x_pct, y, q, ha="right", **kw))
    bracket(a0, a1, INK)
    print(f"list ({G[run[0]]['category']}): " + "; ".join(f"{n}  {q}" for n, q in zip(rows, pcts)))

direct = place_direct(in_list, list(reserved))
for n in order:
    if n in direct:
        (x, y), t, kw, rc, a_lab = direct[n]
        ax.text(x, y, t, **kw)
        texts.append((rc, t))
        if a_lab is not None:
            top = Ri + math.ceil(dots[n]["all"] / arcs[n]["cols"]) * p      # outer edge of its dots
            q0, q1 = xy(arcs[n]["mid"], top + 0.3), xy(a_lab, LAB_R + LEADER_OUT - 0.35)
            ax.plot([q0[0], q1[0]], [q0[1], q1[1]], color=INK, lw=0.35, solid_capstyle="butt")
            print(f"leader: {t}")

labelled = [n for n in order if n in direct or n in prok_names or (n in in_list and G[n]["has_tree"])]
unlabelled = [n for n in order if n not in labelled]

# a run of MIN_RUN or more consecutive groups with no tree and no label of
# their own: one bracket on their arcs and their number, in the grey of a
# group with no tree
MIN_RUN = 3
placed = list(reserved) + [direct[n][3] for n in direct]
run = []
for n in order + [None]:
    if n is not None and n in unlabelled and not G[n]["has_tree"] and not G[n]["is_prokaryotic"] \
            and n not in in_list:
        run.append(n)
        continue
    if len(run) >= MIN_RUN:
        a0, a1 = arcs[run[0]]["a0"], arcs[run[-1]]["a1"]
        t = f"{len(run)} groups, no tree"
        hit = arc_label(t, SOFT, a0, a1, placed, r0=BRACKET_R + 0.9)
        if hit:
            (x, y), t, kw, rc, _ = hit
            ax.text(x, y, t, **kw)
            texts.append((rc, t))
            placed.append(rc)
            bracket(a0, a1, SOFT)
            print(f"run: {t}: {', '.join(run)}")
    run = []

# ── Clade names, at their nodes ────────────────────────────────────────────
# Each name is set at a clear place (no tree line within CLADE_PAD of its
# measured extent, inside the ring, away from the other names): the one
# nearest its node, counting each mm the name's centre lies further in than the
# node, or outside the clade's sector of the ring, as one mm of distance.
# A name more than CLADE_LEADER from its node is tied by a short line to the
# nearer of its node and the arc along its groups.
NODE_R = 0.6
R_CLADE = Ri - 0.85
CLADE_PAD = 0.42
CLADE_LEADER = 3.0
RES = 0.2
# bold: matplotlib reads one face of Helvetica.ttc, so the bold comes from the next faces of the chain
CLADE_KW = dict(fontsize=LABEL_PT, color=INK, fontweight="bold", ha="left", va="bottom",
                fontfamily=["Arial", "DejaVu Sans"])
gx = np.arange(0, FIG_W, RES)
gy = np.arange(0, FIG_H, RES)
busy = np.hypot(gx[None, :] - cx, gy[:, None] - cy) > Ri - 1.0
kpad = int(round(CLADE_PAD / RES))


def block(x0, y0, x1, y1):
    busy[max(0, int(y0 / RES)):int(y1 / RES) + 1, max(0, int(x0 / RES)):int(x1 / RES) + 1] = True


for qx, qy in tree_pts:
    i, j = int(qy / RES), int(qx / RES)
    busy[max(0, i - kpad):i + kpad + 1, max(0, j - kpad):j + kpad + 1] = True
node_xy = [xy(c[1]["_ang"], c[1]["_rad"]) for c in clades]
for nx, ny in node_xy:
    block(nx - NODE_R - 0.5, ny - NODE_R - 0.5, nx + NODE_R + 0.5, ny + NODE_R + 0.5)

for (name, node, members, src, is_clade), (nx, ny) in zip(clades, node_xy):
    t = f"{name}  {pct_text(share(members))}"
    x0, y0, x1, y1 = offsets(t, **CLADE_KW)
    w, h = int(math.ceil((x1 - x0) / RES)), int(math.ceil((y1 - y0) / RES))
    sat = np.zeros((busy.shape[0] + 1, busy.shape[1] + 1), dtype=np.int32)
    sat[1:, 1:] = busy.cumsum(0).cumsum(1)
    full = sat[h:, w:] - sat[:-h, w:] - sat[h:, :-w] + sat[:-h, :-w]     # busy cells under each placement
    ii, jj = np.nonzero(full == 0)
    if not len(ii):
        raise SystemExit(f"no clear place for the clade label {t!r}")
    rx0, ry0 = jj * RES, ii * RES
    dist = np.hypot(np.maximum.reduce([rx0 - nx, np.zeros_like(rx0), nx - rx0 - w * RES]),
                    np.maximum.reduce([ry0 - ny, np.zeros_like(ry0), ny - ry0 - h * RES]))
    mx, my = rx0 + (x1 - x0) / 2 - cx, cy - (ry0 + (y1 - y0) / 2)       # centre of the name, canvas frame
    mr = np.hypot(mx, my)
    a_lo, a_hi = arcs[members[0]]["a0"], arcs[members[-1]]["a1"]
    ma = (np.arctan2(my, mx) - a_lo) % (2 * math.pi) + a_lo             # angle, unwrapped from the sector's start
    out = np.minimum(np.maximum(ma - a_hi, 0), np.maximum(a_lo + 2 * math.pi - ma, 0))
    cost = dist + np.maximum(node["_rad"] - mr, 0) + np.where(ma > a_hi, out * mr, 0)
    k = int(cost.argmin())
    rc = (rx0[k], ry0[k], rx0[k] + (x1 - x0), ry0[k] + (y1 - y0))
    block(rc[0] - 1.0, rc[1] - 1.0, rc[2] + 1.0, rc[3] + 1.0)
    # the clade's extent on the ring: a thin arc between the tips and the dots
    a0, a1 = arcs[members[0]]["a0"] + 0.15 * dth, arcs[members[-1]]["a1"] - 0.15 * dth
    kk = max(2, int((a1 - a0) / 0.01))
    pts = [xy(a0 + (a1 - a0) * i / kk, R_CLADE) for i in range(kk + 1)]
    ax.plot([q[0] for q in pts], [q[1] for q in pts], color=INK, lw=0.5,
            solid_capstyle="butt", zorder=3)
    if is_clade:
        ax.add_patch(Circle((nx, ny), NODE_R, facecolor=INK, edgecolor="white", linewidth=0.3, zorder=4))
    # a name set away from its node is tied to the nearer of its node and its arc
    ends = [(dist[k], nx, ny)] if is_clade else []
    for qx, qy in pts:
        ex, ey = min(max(qx, rc[0]), rc[2]), min(max(qy, rc[1]), rc[3])
        ends.append((math.hypot(qx - ex, qy - ey), qx, qy))
    gap, qx, qy = min(ends)
    if dist[k] > CLADE_LEADER and gap > 1.0:
        ex, ey = min(max(qx, rc[0]), rc[2]), min(max(qy, rc[1]), rc[3])
        f = 0.7 / gap
        ax.plot([ex + (qx - ex) * f, qx], [ey + (qy - ey) * f, qy], color=INK, lw=0.35,
                solid_capstyle="butt", zorder=3)
    ax.text(rc[0] - x0, rc[1] - y0, t, zorder=4, **CLADE_KW)
    texts.append((rc, t))
    print(f"clade {t} ({src}; label {dist[k]:.1f} mm from its node): {', '.join(members)}")

# ── No two texts overlap; nothing leaves the canvas ────────────────────────
for i, (ra, sa) in enumerate(texts):
    if ra[0] < 0 or ra[1] < 0 or ra[2] > FIG_W or ra[3] > FIG_H:
        raise SystemExit(f"text leaves the canvas: {sa!r}")
    for rb, sb in texts[i + 1:]:
        if overlap(ra, rb):
            raise SystemExit(f"texts overlap: {sa!r} / {sb!r}")

afc.savefig(plt, afc.FIGURES_DIR, "figure_tree_of_groups", facecolor="white", pad_inches=0.03)

print(f"SCALE={SCALE} K={K} pitch={p:.3f}mm dot diameter={2 * dotR:.3f}mm eukaryotic dots={n_euk_dots}; "
      f"text {LABEL_PT} pt, headline {head_pt:.2f} pt")
print(f"C0={c0:,.0f} denom={denom:,} pct={100 * c0 / denom:.2f}")
print(f"groups: {n_euk} eukaryotic, {n_tree} with tree, {n_none} without; "
      f"described without tree {desc_none:,} ({100 * desc_none / denom:.1f}%)")
print(f"vert+seed plants named {core_pct:.2f}%  rest {rest_pct:.2f}%  prok {prok_desc:,}")
print("bins:", ", ".join(BINS))
print("texts:", " | ".join(s for _, s in texts))
print(f"labelled {len(labelled)}; unlabelled, clockwise from the top:", ", ".join(unlabelled))
