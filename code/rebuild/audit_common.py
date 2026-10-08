"""Shared helpers for the audit experiments. Everything reads the deposited
atlas files directly; nothing is taken on trust from the manuscript."""
import csv, os, re

ATLAS = os.path.expanduser("~/System/Research/Diversity/phylo-atlas")
FIGS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "figs")

INK, GREY = "#22303C", "#8B96A0"
ACC, GOOD, WARN = "#9E2B25", "#5D8552", "#A07D3E"
BLUE, PURPLE = "#46688F", "#8A76A3"


def p(*a):
    return os.path.join(ATLAS, *a)


def out(name):
    return os.path.abspath(os.path.join(FIGS, name))


def rows(path):
    with open(p(path), newline="") as fh:
        return list(csv.DictReader(fh))


def num(x, default=None):
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return default


def style(plt):
    plt.rcParams.update({
        "font.size": 9, "axes.labelsize": 9.5, "axes.titlesize": 10,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "axes.edgecolor": GREY,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": GREY, "ytick.color": GREY,
    })


# ---------------------------------------------------------------- newick ----
def parse_newick(text):
    """Recursive-descent Newick parser -> nested (children, name, brlen).

    WHITESPACE AT STRUCTURAL POSITIONS. Newick permits whitespace between
    structural tokens, and some deposited trees use it (cetaceans.nwk puts a
    space after every comma). An earlier version of this parser did not skip
    it, so the space became part of the next token's "name" and the tree read
    as 4 leaves instead of 85. Structural whitespace is now removed before
    parsing. Whitespace INSIDE a quoted label is preserved, because a quoted
    label may legitimately contain spaces.
    """
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    # drop whitespace outside quoted labels
    out, quote = [], None
    for ch in s:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            out.append(ch)
        elif ch in " \t\r\n":
            continue
        else:
            out.append(ch)
    s = "".join(out)
    pos = 0

    def node():
        nonlocal pos
        kids = []
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                break
            if pos < len(s) and s[pos] == ")":
                pos += 1
        start = pos
        while pos < len(s) and s[pos] not in "(),":
            pos += 1
        chunk = s[start:pos]
        name, _, bl = chunk.partition(":")
        b = bl.strip()
        if b.lower() in ("nan", "-nan", ""):
            v = float("nan") if b else 0.0
        else:
            v = num(b, 0.0)
            if v is None:
                v = 0.0
        return (kids, name.strip(), v)

    return node()


def leaf_depths(path):
    """Root-to-tip distances for every leaf.

    Tests the paper's own definition of a chronogram (L43): a time-calibrated
    tree is ultrametric, so every leaf sits the same distance from the root.
    Returns (depths, n_leaves, n_distinct_labels, n_nan_branch_lengths).
    """
    root = parse_newick(open(path).read())
    depths, labels, nan_bl = [], [], 0

    stack = [(root, 0.0)]
    while stack:
        (kids, name, bl), acc = stack.pop()
        import math
        if bl != bl:                      # NaN branch length
            nan_bl += 1
            bl = 0.0
        d = acc + bl
        if kids:
            for k in kids:
                stack.append((k, d))
        else:
            depths.append(d)
            labels.append(name)
    return depths, len(depths), len(set(labels)), nan_bl
