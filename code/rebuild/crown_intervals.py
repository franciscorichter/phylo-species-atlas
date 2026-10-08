#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""crown_intervals.py -- crown-age intervals of the canonical trees, read from
what the source studies deposited (finding S02; decisions D5, D7, D17).

The crown age of a partition is the root height of its canonical file: after
decision D7 the root of every canonical tree is the most recent common
ancestor of the tips retained in the partition.  rebuild_atlas.py measures that
point from the shipped file.  This script measures the INTERVAL around it from
the deposit named in

    Data/crown_interval_sources.csv   one row per deposit.  Columns:
        partition      the partition whose canonical tree the row describes
        source         path(s) under Data/ ('|'-separated; a '*' glob is
                       expanded and sorted).  'standardized/trees/<file>'
                       reads the shipped canonical file itself.
        kind           what the deposit is (see KINDS)
        outside        source labels outside the focal clade, ';'-separated:
                       the tips decision D7 prunes, plus any listed in
                       focal_outgroups.  The focal node in every tree of the
                       source is the MRCA of all other labels.
        focal_outgroups  the subset of `outside` that belongs to the
                       partition and stays in the canonical file (an outgroup
                       of the source's ingroup that is still in the partition,
                       e.g. Antipatharia for the Scleractinia study, D21)
        annotation     for node_hpd_annotation, the annotation key holding
                       the interval (e.g. height_95%_HPD) and, after '|', the
                       key holding the deposit's own point summary
        scale          multiplier to Ma (1, or 100 where the source is in
                       units of 100 Myr)
        burnin         fraction of each sample file dropped from its start
                       (non-zero only for raw MCMC logs)
        evidence       why the kind and the parameters are what they are

KINDS, and whether each makes an interval RECOVERABLE (rung C1b):
    node_hpd_annotation       yes  per-node interval written in the summary
                                   tree by the dating software (TreeAnnotator,
                                   MCMCTree)
    node_label_interval       yes  per-node 'upper_lower' internal labels
                                   (PhyloBayes readdiv)
    posterior_sample          yes  dated trees sampled from the posterior
    bootstrap_dating_sample   yes  each bootstrap tree dated with the same
                                   calibrations (treePL)
    calibration_prior_sample  no   the source's tree set dated under root
                                   ages drawn from the calibration density;
                                   the spread is the prior, not an estimate,
                                   so the interval is set by a calibration
                                   (c1b_exclusion prior_sampled, D29, D51)
    taxon_addition_sample     no   one dated backbone, with taxonomically
                                   placed tips added at random per tree; the
                                   spread is placement, not dating
    date_assignment_sample    no   per-node published dates drawn at random
                                   onto one topology; the spread is between
                                   literature dates, not an estimate

For a sample the interval is the equal-tailed 95% interval of the focal node's
height across trees (2.5 and 97.5% quantiles, linear interpolation), with the
median as its centre.  It is not a highest-posterior-density interval.  For an
annotated or labelled summary tree it is the interval written on the focal
node.  Every row also records the focal node's branch-length height in the
first (or only) tree of the source, which verify_atlas.py compares with the
canonical file's root height: that equality is what shows the focal node is
the canonical crown.

DEGENERATE.  An interval narrower than DEGENERATE_WIDTH_MYR is flagged
degenerate (decision D5): the node is fixed, typically by a calibration bound.
AT BOUND.  A centre (median) within BOUND_TOL_MYR of the interval's lower or
upper end means at least 47.5% of the mass sits on one value: the node is held
at a hard calibration bound (Diatoms: 190 Ma, the diatomMRCA minimum).
Both flags are computed for every row by the same rule; no row is dropped.
PRIOR REPRODUCED (decision D25).  Data/calibration_register.csv records,
for every row whose kind counts, the calibration the source placed on the
focal node.  Where it is two-sided (a lower and an upper bound, or a
parametric density with both 95% quantiles set by the calibration), its 95%
prior range [p_lo, p_hi] is compared with the reported interval [lo, hi]:
    containment = |[lo, hi] & [p_lo, p_hi]| / (hi - lo)
    width_ratio = (hi - lo) / (p_hi - p_lo)
The interval reproduces the prior when containment >= PRIOR_CONTAINMENT_MIN
(the interval lies inside the prior's 95% range) and width_ratio >=
PRIOR_WIDTH_RATIO_MIN (the data did not even halve its width).  The prior
range is the effective (marginal) prior from a run of the source's own input
file sampling from the prior where one exists (prior_basis
effective_prior_run), otherwise the calibration as specified.  A one-sided,
absent or unchecked calibration is recorded and does not exclude by itself.
BAND (round F4, F4-09; decision D52).  Where the crown carries no two-sided
calibration of its own, calibrations on the nodes around it can still
confine it: the crown is no younger than the oldest minimum on it or on a
node below it, and no older than the youngest maximum on it or on a node
above it.  The register records that band (band_lo_ma, band_hi_ma, with
child_calibration and band_evidence), and the same test is applied to it
(prior_range = band); with no band the interval is compared with no prior.
A row counts for C1b only when its kind counts AND it is neither degenerate
nor at a bound (decision D22) nor reproduces its prior (D25); c1b_exclusion
names the reason otherwise.

OUTPUT: Data/crown_intervals.csv, one row per register row.  Rows whose source
is not present locally keep their committed values; --check reports them.

USAGE
    python3 scripts/crown_intervals.py [--only PARTITION ...] [--jobs N]
                                       [--data-root DIR] [--check]
Standard library only, so it runs on any compute host.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
REGISTER = os.path.join(ROOT, "Data", "crown_interval_sources.csv")
CALIBRATIONS = os.path.join(ROOT, "Data", "calibration_register.csv")
OUT = os.path.join(ROOT, "Data", "crown_intervals.csv")

KINDS = {
    "node_hpd_annotation": True,
    "node_label_interval": True,
    "posterior_sample": True,
    "bootstrap_dating_sample": True,
    "calibration_prior_sample": False,
    "taxon_addition_sample": False,
    "date_assignment_sample": False,
}
PRIOR_SAMPLED_KINDS = ("calibration_prior_sample",)   # D29
SAMPLE_KINDS = ("posterior_sample", "bootstrap_dating_sample",
                "calibration_prior_sample", "taxon_addition_sample",
                "date_assignment_sample")
QUANTILES = (0.025, 0.5, 0.975)        # equal-tailed 95% interval and median
DEGENERATE_WIDTH_MYR = 1e-3            # D5 rule parameter, applied to every row
BOUND_TOL_MYR = 1e-2                   # D22 rule parameter, applied to every row
PRIOR_CONTAINMENT_MIN = 0.9            # D25 rule parameters, applied to every
PRIOR_WIDTH_RATIO_MIN = 0.5            # row with a two-sided calibration
CAL_CLASSES = ("two_sided", "one_sided", "none", "not_checked")

OUT_HEADER = ["partition", "kind", "counts_for_c1b", "source",
              "source_sha256", "n_files", "n_trees", "burnin_trees",
              "n_focal_labels", "n_outside_labels", "source_point_height_ma",
              "interval_lo_ma", "interval_center_ma", "interval_hi_ma",
              "center_statistic", "interval_width_myr", "degenerate",
              "at_bound", "c1b_exclusion", "n_focal_outgroups",
              "source_root_height_ma", "calibration", "prior_lo_ma",
              "prior_hi_ma", "prior_containment", "prior_width_ratio",
              "prior_reproduced", "prior_q_lo", "prior_q_center",
              "prior_q_hi", "prior_range"]

csv.field_size_limit(1 << 30)
TOKEN_RE = re.compile(r"\[[^\]]*\]|'(?:[^']|'')*'|[(),;:]|[^(),;:\['\s]+")
COMMENT_RE = re.compile(r"\[[^\]]*\]")


# --------------------------------------------------------------------------
# reading trees
# --------------------------------------------------------------------------

def iter_tree_strings(path):
    """Every tree statement of a Newick or NEXUS file, translate table
    applied to the leaf labels (by the parser), comments kept."""
    translate, in_tr, nexus, buf = {}, False, None, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if nexus is None and s:
                nexus = s.upper().startswith("#NEXUS")
            if buf is not None:
                buf += s
                if s.endswith(";"):
                    yield buf, translate
                    buf = None
                continue
            if not s:
                continue
            if nexus:
                low = s.lower()
                if low.startswith("translate"):
                    in_tr = True
                    s = s[len("translate"):].strip()
                    if not s:
                        continue
                if in_tr:
                    for item in s.rstrip(";").split(","):
                        kv = item.split(None, 1)
                        if len(kv) == 2:
                            translate[kv[0]] = kv[1].strip().strip("'")
                    if s.endswith(";"):
                        in_tr = False
                    continue
                # a bracketed comment between the name and '=' (BEAST:
                # 'tree STATE_1 [&posterior=...] = [&R] (...)') may hold '='
                m = re.match(r"^u?tree\b(?:[^=\[]|\[[^\]]*\])*=\s*"
                             r"(?:\[[^\]]*\]\s*)*(\(.*)$", s, re.I)
                if m:
                    buf = m.group(1)
                    if buf.endswith(";"):
                        yield buf, translate
                        buf = None
            elif s.startswith("("):
                buf = s
                if buf.endswith(";"):
                    yield buf, translate
                    buf = None


class Tree(object):
    """Nodes in preorder: parent[i] < i for every non-root node."""
    __slots__ = ("parent", "bl", "label", "comment", "leaf")


def parse(text, translate, keep_comments):
    if not keep_comments:
        text = COMMENT_RE.sub("", text)
    t = Tree()
    parent, bl, label, comment, leaf = [-1], [0.0], [""], [""], [False]
    cur = 0
    expect_bl = False
    for m in TOKEN_RE.finditer(text):
        tok = m.group(0)
        if tok[0] == "[":
            if keep_comments:
                comment[cur] += tok
            continue
        if tok == "(":
            parent.append(cur); bl.append(0.0); label.append("")
            comment.append(""); leaf.append(False)
            cur = len(parent) - 1
            expect_bl = False
        elif tok == ",":
            p = parent[cur]
            parent.append(p); bl.append(0.0); label.append("")
            comment.append(""); leaf.append(False)
            cur = len(parent) - 1
            expect_bl = False
        elif tok == ")":
            cur = parent[cur]
            expect_bl = False
        elif tok == ":":
            expect_bl = True
        elif tok == ";":
            break
        elif expect_bl:
            try:
                bl[cur] += float(tok)
            except ValueError:
                pass
            expect_bl = False
        else:
            if tok.startswith("'"):
                tok = tok[1:-1].replace("''", "'")
            label[cur] = tok
    # a node with no child is a leaf
    haskid = [False] * len(parent)
    for i in range(1, len(parent)):
        haskid[parent[i]] = True
    for i in range(len(parent)):
        leaf[i] = not haskid[i]
        if leaf[i] and translate:
            label[i] = translate.get(label[i], label[i])
    t.parent, t.bl, t.label, t.comment, t.leaf = parent, bl, label, comment, leaf
    return t


def focal_node(t, outside):
    """(node index, its branch-length height, number of focal leaves).
    The focal node is the MRCA of every leaf whose label is not in `outside`;
    its height is the deepest leaf below it minus its own depth."""
    n = len(t.parent)
    depth = [0.0] * n
    for i in range(1, n):
        depth[i] = depth[t.parent[i]] + t.bl[i]
    cnt = [0] * n
    deep = [-1.0] * n
    total = 0
    # (the root height of the whole tree is deep[0] - depth[0])
    for i in range(n):
        if t.leaf[i]:
            deep[i] = depth[i]
            if t.label[i] not in outside:
                cnt[i] = 1
                total += 1
    for i in range(n - 1, 0, -1):
        p = t.parent[i]
        cnt[p] += cnt[i]
        if deep[i] > deep[p]:
            deep[p] = deep[i]
    best = 0
    for i in range(n):              # preorder: the last full node is deepest
        if cnt[i] == total and depth[i] >= depth[best]:
            best = i
    return best, deep[best] - depth[best], total


def root_height(t):
    """Deepest leaf depth below the root."""
    n = len(t.parent)
    depth = [0.0] * n
    for i in range(1, n):
        depth[i] = depth[t.parent[i]] + t.bl[i]
    return max(depth[i] for i in range(n) if t.leaf[i])


def parse_annotation(comment, key):
    """Value(s) of `key` in a BEAST/MCMCTree comment: a float or a
    (lo, hi) pair."""
    body = comment.replace("[&", ",").replace("]", "")
    m = re.search(r"(?:^|,)\s*" + re.escape(key) + r"\s*=\s*(\{[^}]*\}|[^,]+)",
                  body)
    if not m:
        return None
    v = m.group(1).strip()
    if v.startswith("{"):
        a = [float(x) for x in v.strip("{}").split(",")]
        return (min(a), max(a))
    return float(v)


def quantile(xs, q):
    s = sorted(xs)
    if not s:
        return None
    pos = (len(s) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def sha256_of(paths):
    """sha256 of one file, or of the sorted per-file digests for several."""
    digs = []
    for p in paths:
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        digs.append(h.hexdigest())
    if len(digs) == 1:
        return digs[0]
    return "set:" + hashlib.sha256("".join(sorted(digs)).encode()).hexdigest()


def _sample_file(args):
    path, outside, burnin, scale = args
    trees = list(iter_tree_strings(path))
    drop = int(round(burnin * len(trees)))
    heights, first, nfocal = [], None, None
    for k, (s, tr) in enumerate(trees):
        if k < drop:
            continue
        t = parse(s, tr, False)
        _i, h, nf = focal_node(t, outside)
        heights.append(h * scale)
        if first is None:
            first, nfocal = h * scale, nf
            root_h = root_height(t) * scale
    return heights, first, nfocal, len(trees), drop, root_h


# --------------------------------------------------------------------------
# one register row
# --------------------------------------------------------------------------

def measure(row, data_root, jobs):
    kind = row["kind"]
    if kind not in KINDS:
        sys.exit("unknown kind %r for %s" % (kind, row["partition"]))
    paths = []
    for spec in row["source"].split("|"):
        spec = spec.strip()
        hit = sorted(glob.glob(os.path.join(data_root, spec)))
        if not hit:
            return None
        paths += hit
    outside = {x.strip() for x in row["outside"].split(";") if x.strip()}
    kept = {x.strip() for x in (row.get("focal_outgroups") or "").split(";")
            if x.strip()}
    if not kept <= outside:
        sys.exit("%s: focal_outgroups must be a subset of outside"
                 % row["partition"])
    scale = float(row["scale"] or 1)
    burnin = float(row["burnin"] or 0)
    out = dict(partition=row["partition"], kind=kind,
               counts_for_c1b="yes" if KINDS[kind] else "no",
               source=row["source"], source_sha256=sha256_of(paths),
               n_files=len(paths), n_outside_labels=len(outside),
               n_focal_outgroups=len(kept))
    if kind in SAMPLE_KINDS:
        args = [(p, outside, burnin, scale) for p in paths]
        if jobs > 1 and len(args) > 1:
            with ProcessPoolExecutor(min(jobs, len(args))) as ex:
                res = list(ex.map(_sample_file, args))
        else:
            res = [_sample_file(a) for a in args]
        heights = [h for r in res for h in r[0]]
        lo, med, hi = (quantile(heights, q) for q in QUANTILES)
        out.update(n_trees=sum(r[3] for r in res) - sum(r[4] for r in res),
                   burnin_trees=sum(r[4] for r in res),
                   n_focal_labels=res[0][2], source_point_height_ma=res[0][1],
                   interval_lo_ma=lo, interval_center_ma=med,
                   interval_hi_ma=hi, center_statistic="median",
                   source_root_height_ma=res[0][5])
    else:
        s, tr = next(iter_tree_strings(paths[0]))
        keep = kind == "node_hpd_annotation"
        t = parse(s, tr, keep)
        i, h, nf = focal_node(t, outside)
        if kind == "node_hpd_annotation":
            ikey, _, ckey = row["annotation"].partition("|")
            iv = parse_annotation(t.comment[i], ikey.strip())
            cv = parse_annotation(t.comment[i], ckey.strip()) \
                if ckey.strip() else None
            if not isinstance(iv, tuple):
                sys.exit("%s: no %s on the focal node" % (row["partition"],
                                                          ikey))
            lo, hi = iv[0] * scale, iv[1] * scale
            center = cv * scale if isinstance(cv, float) else ""
            stat = ckey.strip()
        else:
            m = re.match(r"^\s*([0-9.eE+-]+)_([0-9.eE+-]+)\s*$", t.label[i])
            if not m:
                sys.exit("%s: focal node label %r is not upper_lower"
                         % (row["partition"], t.label[i]))
            a, b = float(m.group(1)) * scale, float(m.group(2)) * scale
            lo, hi, center, stat = min(a, b), max(a, b), "", ""
        out.update(n_trees=1, burnin_trees=0, n_focal_labels=nf,
                   source_point_height_ma=h * scale, interval_lo_ma=lo,
                   interval_center_ma=center, interval_hi_ma=hi,
                   center_statistic=stat,
                   source_root_height_ma=root_height(t) * scale)
    out.update(flags(kind, out["interval_lo_ma"], out["interval_center_ma"],
                     out["interval_hi_ma"], out["source_point_height_ma"],
                     load_calibrations().get(row["partition"])))
    return out


def prior_range(cal):
    """D25, D52.  (p_lo, p_hi, source) of the range rule (p) compares an
    interval with: the crown's own two-sided calibration ('calibration'),
    otherwise the band the calibrations on the nodes around it allow
    ('band'), otherwise none."""
    if not cal:
        return None, None, ""
    if cal["calibration"] == "two_sided":
        return float(cal["prior_lo_ma"]), float(cal["prior_hi_ma"]), "calibration"
    blo = (cal.get("band_lo_ma") or "").strip()
    bhi = (cal.get("band_hi_ma") or "").strip()
    if blo and bhi:
        if float(blo) >= float(bhi):
            sys.exit("calibration band %s-%s is not ordered" % (blo, bhi))
        return float(blo), float(bhi), "band"
    return None, None, ""


def prior_test(lo, hi, cal):
    """D25.  (calibration, p_lo, p_hi, containment, width_ratio, reproduced)
    for one interval and its calibration-register row (None: no row)."""
    if not cal:
        return "", None, None, None, None, False
    klass = cal["calibration"]
    if klass not in CAL_CLASSES:
        sys.exit("unknown calibration class %r" % klass)
    plo, phi, _src = prior_range(cal)
    if plo is None:
        return klass, None, None, None, None, False
    w = hi - lo
    ov = max(0.0, min(hi, phi) - max(lo, plo))
    cont = ov / w if w > 0 else 1.0
    ratio = w / (phi - plo)
    return (klass, plo, phi, cont, ratio,
            cont >= PRIOR_CONTAINMENT_MIN and ratio >= PRIOR_WIDTH_RATIO_MIN)


PRIOR_DENSITY_RE = re.compile(r"^\s*(\w+)\(([^)]*)\)\s*$")
SOFT_TAIL = 0.025                      # Yang & Rannala 2006 soft bounds: 2.5%
                                       # of the prior outside each bound


def prior_quantile(spec, x):
    """Round F3 (F3-14).  The quantile of x under the calibration density
    the source specifies on the focal node (register column prior_density):
        uniform(lo,hi)            uniform between hard bounds
        soft_uniform(lo,hi)       uniform inside soft bounds holding 95% of
                                  the prior (PhyloBayes -sb); outside the
                                  bounds only the side is reported
        offset_exp_q95(min,max)   exponential offset to min whose 95%
                                  quantile is max (Upham et al. 2019 NDexp)
        offset_gamma(off,a,s)     gamma with shape a and scale s offset by off
    Returns a string: the quantile to three decimals, '<0.025' or '>0.975'
    beyond soft bounds, '' without a density."""
    if not spec:
        return ""
    m = PRIOR_DENSITY_RE.match(spec)
    if not m:
        sys.exit("unreadable prior_density %r" % spec)
    form, args = m.group(1), [float(v) for v in m.group(2).split(",")]
    if form == "uniform":
        lo, hi = args
        q = min(1.0, max(0.0, (x - lo) / (hi - lo)))
    elif form == "soft_uniform":
        lo, hi = args
        if x < lo:
            return "<%.3f" % SOFT_TAIL
        if x > hi:
            return ">%.3f" % (1 - SOFT_TAIL)
        q = SOFT_TAIL + (1 - 2 * SOFT_TAIL) * (x - lo) / (hi - lo)
    elif form == "offset_exp_q95":
        lo, hi = args
        rate = -math.log(0.05) / (hi - lo)
        q = 0.0 if x <= lo else 1.0 - math.exp(-rate * (x - lo))
    elif form == "offset_gamma":
        off, a, sc = args
        if a != int(a):
            sys.exit("offset_gamma needs an integer shape, got %r" % spec)
        z = max(0.0, (x - off) / sc)
        # regularised lower incomplete gamma for integer shape (Erlang)
        q = 1.0 - math.exp(-z) * sum(z ** k / math.factorial(k)
                                     for k in range(int(a)))
    else:
        sys.exit("unknown prior_density form %r" % form)
    return "%.3f" % q


def prior_quantiles(lo, center, hi, point, cal):
    spec = (cal or {}).get("prior_density", "") or ""
    c = center if isinstance(center, float) else point
    return dict(prior_q_lo=prior_quantile(spec, lo),
                prior_q_center=prior_quantile(spec, c)
                if c is not None else "",
                prior_q_hi=prior_quantile(spec, hi))


def load_calibrations(path=CALIBRATIONS):
    if not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["partition"]: r for r in csv.DictReader(fh)}


def flags(kind, lo, center, hi, point, cal=None):
    """width, degenerate, at_bound, the D25 prior test, c1b_exclusion and
    counts_for_c1b for one interval, by the rules above (D5, D22, D25).
    Shared with verify_atlas.py."""
    w = hi - lo
    c = center if isinstance(center, float) else point
    degenerate = w < DEGENERATE_WIDTH_MYR
    at_bound = (not degenerate and c is not None
                and (c - lo < BOUND_TOL_MYR or hi - c < BOUND_TOL_MYR))
    klass, plo, phi, cont, ratio, repro = prior_test(lo, hi, cal)
    # D29, D51: a tree set dated under root ages drawn from the calibration
    # prior returns that prior at the crown by construction; it is set by a
    # calibration, under its own rule (prior_sampled), not by the D25 test.
    why = ("prior_sampled" if kind in PRIOR_SAMPLED_KINDS
           else "kind" if not KINDS[kind] else "degenerate" if degenerate
           else "at_bound" if at_bound
           else "prior_reproduced" if repro else "")
    return dict(interval_width_myr=w,
                degenerate="yes" if degenerate else "no",
                at_bound="yes" if at_bound else "no",
                calibration=klass, prior_lo_ma=plo, prior_hi_ma=phi,
                prior_containment=cont, prior_width_ratio=ratio,
                prior_reproduced="yes" if repro else "no",
                **prior_quantiles(lo, center, hi, point, cal),
                prior_range=prior_range(cal)[2],
                c1b_exclusion=why,
                counts_for_c1b="yes" if not why else "no")


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def fmt(v):
    if isinstance(v, float):
        return "%.6f" % v
    return "" if v is None else str(v)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--data-root", default=os.path.join(ROOT, "Data"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--prior-log", metavar="BEAST_LOG",
                    help="print the equal-tailed 95%% interval and median of "
                         "--column in a BEAST trace log (a sampleFromPrior "
                         "run) after --log-burnin, by the same quantile rule; "
                         "the effective prior of Data/calibration_register.csv")
    ap.add_argument("--column", default="Tree.height")
    ap.add_argument("--log-burnin", type=float, default=0.1)
    a = ap.parse_args(argv)
    if a.prior_log:
        with open(a.prior_log, encoding="utf-8") as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh
                    if ln.strip() and not ln.startswith("#")]
        head, body = rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]
        i = head.index(a.column)
        burn = int(len(body) * a.log_burnin)
        xs = [float(r[i]) for r in body[burn:]]
        print("samples=%d burnin=%d kept=%d last_state=%s  2.5%%=%.4f  "
              "50%%=%.4f  97.5%%=%.4f" % (len(body), burn, len(xs),
                                          body[-1][0], *(quantile(xs, q)
                                                         for q in QUANTILES)))
        return 0
    with open(REGISTER, newline="", encoding="utf-8") as fh:
        reg = list(csv.DictReader(fh))
    old = {}
    if os.path.exists(a.out):
        with open(a.out, newline="", encoding="utf-8") as fh:
            old = {r["partition"] + "|" + r["kind"]: r
                   for r in csv.DictReader(fh)}
    rows, bad = [], 0
    for r in reg:
        key = r["partition"] + "|" + r["kind"]
        if a.only and r["partition"] not in a.only:
            if key in old:
                rows.append(old[key])
            continue
        got = measure(r, a.data_root, a.jobs)
        if got is None:
            print("%-13s %-25s source not present here; committed row kept"
                  % (r["partition"], r["kind"]))
            if key in old:
                rows.append(old[key])
            else:
                bad += 1
            continue
        new = {k: fmt(got.get(k, "")) for k in OUT_HEADER}
        if a.check and key in old:
            diff = [k for k in OUT_HEADER if old[key].get(k, "") != new[k]]
            if diff:
                bad += 1
                print("%-13s CHANGED: %s" % (r["partition"], diff))
        print("%-13s %-25s n=%s focal=%s point=%s interval=[%s, %s]"
              % (r["partition"], r["kind"], new["n_trees"],
                 new["n_focal_labels"], new["source_point_height_ma"],
                 new["interval_lo_ma"], new["interval_hi_ma"]))
        rows.append(new)
    # D25 and the D5/D22 flags are re-applied to EVERY row, including rows
    # kept from the committed file because their source is not present here:
    # they are functions of the committed interval and the register.
    cals = load_calibrations()
    for row in rows:
        f = flags(row["kind"], float(row["interval_lo_ma"]),
                  _num(row["interval_center_ma"]), float(row["interval_hi_ma"]),
                  _num(row["source_point_height_ma"]), cals.get(row["partition"]))
        row.update({k: fmt(v) for k, v in f.items()})
    if a.check:
        print("%d row(s) differ or cannot be produced" % bad)
        return 1 if bad else 0
    rows.sort(key=lambda r: (r["partition"], r["kind"]))
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, OUT_HEADER, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("wrote %s (%d rows)" % (os.path.relpath(a.out, ROOT), len(rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
