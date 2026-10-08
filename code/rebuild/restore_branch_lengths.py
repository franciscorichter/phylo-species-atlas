#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""restore_branch_lengths.py -- undo the digit-remap corruption of the shipped trees.

WHAT HAPPENED
    Between the 2026-04-09 checkpoint and the 2026-05-12 submission copy, the
    shipped trees were relabelled from the 641,763-row dictionary to the
    637,619-row dictionary by a step that rewrote digit runs that were not tip
    labels.  The rule, which reproduces all 257 archived corrupted copies byte
    for byte from their pre-remap copies, is: every run of digits matching
    \b\d+ that is immediately followed by ',', ')' or ':' had its leading zeros
    dropped and, when the resulting integer was an id of the old dictionary,
    was replaced by the new id of the same standardised name.  Runs followed
    by anything else were left untouched: root-edge lengths (followed by ';'),
    mantissas before an exponent, quoted labels and internal labels that begin
    with a letter.  Affected branch lengths and internal labels:

        0.007143      -> 0.7146        (7143 is an old id; new id 7146)
        27.05691095   -> 27.5691095    (no old id 5691095; leading zero dropped)
        9.0432274     -> 9.435619      (old id 432274 is new id 435619)

    The old -> new tip map matches names after replacing spaces by
    underscores; for bryozoa (720 tips) and frogs (5,326 tips) the names
    differ only in that respect, so the match is not an exact string match
    there.

    The pre-remap copies survive at
    Paper/archived/wrong_root_standardized_20260512_012852/trees/ (the
    "wrong_root" in the name refers to the repository root directory they were
    copied from, not to tree rooting).  For every shipped tree the pre-remap
    copy has the identical token structure, and every tip label maps exactly
    through old id -> standardised name -> new id.

WHAT THIS SCRIPT DOES
    For every tree in Data/standardized/trees/ (and every site copy that is
    byte-identical to a shipped tree), it rebuilds the file from the pre-remap
    copy, replacing ONLY the tip labels (old id -> new id).  Every other token --
    branch lengths, internal node labels, comments -- is copied verbatim from
    the pre-remap file.  It refuses to write a tree unless

      * the two files have the same token structure, and
      * every tip of the pre-remap file maps to exactly the tip the shipped file
        carries at the same position.

    The file being replaced is MOVED to
    Data/_archive/2026-10-revision/corrupted_branch_lengths/<same relative path>
    (never deleted).  Idempotent: a tree that already equals its restoration is
    left alone.

    A per-tree record is written to
    Data/_archive/2026-10-revision/corrupted_branch_lengths/RESTORED.csv.

USAGE
    python3 scripts/restore_branch_lengths.py [--dry-run]
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
PRE = os.path.join(ROOT, "Paper", "archived",
                   "wrong_root_standardized_20260512_012852")
PRE_TREES = os.path.join(PRE, "trees")
PRE_DICT = os.path.join(PRE, "dictionary.csv")
STD_TREES = os.path.join(ROOT, "Data", "standardized", "trees")
STD_DICT = os.path.join(ROOT, "Data", "standardized", "dictionary.csv")
SITE_PARTS = os.path.join(ROOT, "Data", "site", "data", "partitions")
ARCHIVE = os.path.join(ROOT, "Data", "_archive", "2026-10-revision",
                       "corrupted_branch_lengths")
RECORD = os.path.join(ARCHIVE, "RESTORED.csv")

csv.field_size_limit(1 << 30)

# One Newick token per match: punctuation, ':length', [comment], 'quoted' or bare label.
TOKEN_RE = re.compile(r"[(),;]|:[^,();\[]*|\[[^\]]*\]|'(?:[^']|'')*'|[^(),:;\[']+")


def rel(p):
    return os.path.relpath(p, ROOT)


def load_dict(path):
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            out[row["id"]] = row["standardized_name"]
    return out


def norm_name(s):
    return s.replace(" ", "_")


def tokens(text):
    """[(kind, start, end)] with kind in p (punct), b (branch length),
    c (comment), l (label).  Whitespace-only label tokens are dropped."""
    out = []
    for m in TOKEN_RE.finditer(text):
        t = m.group(0)
        if t in "(),;" and len(t) == 1:
            out.append(("p", m.start(), m.end()))
        elif t.startswith(":"):
            out.append(("b", m.start(), m.end()))
        elif t.startswith("["):
            out.append(("c", m.start(), m.end()))
        elif t.strip():
            out.append(("l", m.start(), m.end()))
    return out


def tip_positions(text, toks):
    """Indices of label tokens that are tips (preceded by '(' or ',')."""
    tips = set()
    prev = None
    for i, (k, s, e) in enumerate(toks):
        if k == "l" and prev in ("(", ",", None):
            tips.add(i)
        if k == "p":
            prev = text[s:e]
        elif k in ("l", "b"):
            prev = k
    return tips


def restore(pre_text, shipped_text, old2new):
    """Return (restored_text, stats) or raise ValueError."""
    tp, ts = tokens(pre_text), tokens(shipped_text)
    if len(tp) != len(ts):
        raise ValueError("token count differs (%d vs %d)" % (len(tp), len(ts)))
    tips = tip_positions(pre_text, tp)
    pieces, last = [], 0
    n_bl = n_bl_diff = n_int_diff = 0
    for i, ((kp, sp, ep), (ks, ss, es)) in enumerate(zip(tp, ts)):
        if kp != ks:
            raise ValueError("token kind differs at token %d" % i)
        a, b = pre_text[sp:ep], shipped_text[ss:es]
        if kp == "p" and a != b:
            raise ValueError("punctuation differs at token %d" % i)
        if kp == "l" and i in tips:
            lab = a.strip().strip("'")
            new = old2new.get(lab)
            if new is None:
                raise ValueError("pre-remap tip %r has no new id" % lab)
            if new != b.strip().strip("'"):
                raise ValueError("pre-remap tip %r maps to %s, shipped has %r"
                                 % (lab, new, b))
            pieces.append(pre_text[last:sp])
            pieces.append(b.strip())
            last = ep
        elif kp == "b":
            n_bl += 1
            if a != b:
                n_bl_diff += 1
        elif kp == "l" and a != b:
            n_int_diff += 1
    pieces.append(pre_text[last:])
    return "".join(pieces), {"branch_lengths": n_bl,
                             "branch_lengths_restored": n_bl_diff,
                             "internal_labels_restored": n_int_diff}


def move_to_archive(path, dry):
    dest = os.path.join(ARCHIVE, rel(path))
    if os.path.exists(dest):
        return dest
    if not dry:
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy2(path, dest)
    return dest


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    old = load_dict(PRE_DICT)
    new = load_dict(STD_DICT)
    name2new = {}
    for k, v in new.items():
        name2new.setdefault(norm_name(v), k)
    old2new = {k: name2new.get(norm_name(v)) for k, v in old.items()}

    shipped = sorted(fn for fn in os.listdir(STD_TREES) if fn.endswith(".nwk"))
    corrupted_text = {}
    records, failures = [], []
    for fn in shipped:
        path = os.path.join(STD_TREES, fn)
        pre_path = os.path.join(PRE_TREES, fn)
        if not os.path.exists(pre_path):
            failures.append((fn, "no pre-remap copy"))
            continue
        with open(path, encoding="utf-8") as fh:
            cur = fh.read()
        with open(pre_path, encoding="utf-8") as fh:
            pre = fh.read()
        try:
            out, st = restore(pre, cur, old2new)
        except ValueError as exc:
            failures.append((fn, str(exc)))
            continue
        changed = out != cur
        if changed:
            corrupted_text[cur] = out
            move_to_archive(path, args.dry_run)
            if not args.dry_run:
                tmp = path + ".tmp-%d" % os.getpid()
                with open(tmp, "w", encoding="utf-8", newline="") as fh:
                    fh.write(out)
                os.replace(tmp, path)
        records.append(dict(file=rel(path), reference=rel(pre_path),
                            changed="yes" if changed else "no", **st))

    # Site copies that are byte-identical to a corrupted shipped tree.
    for dirpath, _dirs, files in os.walk(SITE_PARTS):
        for fn in files:
            if not fn.endswith(".nwk"):
                continue
            path = os.path.join(dirpath, fn)
            with open(path, encoding="utf-8", errors="replace") as fh:
                cur = fh.read()
            if cur in corrupted_text:
                move_to_archive(path, args.dry_run)
                if not args.dry_run:
                    with open(path, "w", encoding="utf-8", newline="") as fh:
                        fh.write(corrupted_text[cur])
                records.append(dict(file=rel(path), reference="(copy of the "
                                    "restored shipped tree)", changed="yes",
                                    branch_lengths="", branch_lengths_restored="",
                                    internal_labels_restored=""))

    if failures:
        for fn, why in failures:
            sys.stderr.write("REFUSED %s: %s\n" % (fn, why))
        sys.exit(1)

    if not args.dry_run and records:
        prior = []
        if os.path.exists(RECORD):
            with open(RECORD, newline="", encoding="utf-8") as fh:
                prior = list(csv.DictReader(fh))
        # Keep the first record of each file: a re-run finds nothing to restore
        # and must not overwrite the record of what was restored.
        by_file = {r["file"]: r for r in prior}
        for r in records:
            if r["file"] not in by_file or by_file[r["file"]]["changed"] == "no":
                by_file[r["file"]] = r
        os.makedirs(ARCHIVE, exist_ok=True)
        hdr = ["file", "reference", "changed", "branch_lengths",
               "branch_lengths_restored", "internal_labels_restored"]
        with open(RECORD, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=hdr, lineterminator="\n")
            w.writeheader()
            for k in sorted(by_file):
                w.writerow({h: by_file[k].get(h, "") for h in hdr})

    n_changed = sum(1 for r in records if r["changed"] == "yes")
    print("%d files examined, %d restored%s" % (len(records), n_changed,
                                                " (dry run)" if args.dry_run else ""))


if __name__ == "__main__":
    main()
