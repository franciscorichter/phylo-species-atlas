#!/usr/bin/env python3
"""The licence of every tree file of the data release, and the files it withholds.

The data release carries tree files in four places: the trees of the catalogue
(``standardized/trees/``), the unpruned source of each reference tree
(``standardized/canonical_sources/``), the former reference trees
(``standardized/superseded_canonicals/``) and the sub-clade trees of the
browser site (``site/data/partitions/*/sub-phylos/``).  ``Data/tree_licences.csv``
records, for each source, the licence its archive states, the evidence and the
redistribution class.  This module joins every candidate file to its row:

    open, open-share-alike   the file is redistributed under the source's licence
    no-licence-stated        the file is withheld: the paper keeps it, the
                             release lists it in registers/withheld_trees.csv
                             with the source to download the original from

A candidate file that matches no row, or a row of another class, stops the
build.  The outputs, all written by ``scripts/rebuild_atlas.py`` and checked
by ``scripts/verify_atlas.py`` (RELEASE):

    Data/withheld_trees.csv          one row per withheld file: its release
                                     path, study, source, licence, the sha256,
                                     tip count, root height and root-to-tip CV
                                     of the file the paper used, and how that
                                     file was derived from the source
    Data/withheld_tips/<release path>.txt
                                     the tip labels of each withheld file,
                                     sorted (a list of species, not a tree);
                                     the release ships them gzipped under
                                     registers/withheld_tips/
    Data/LICENSE-trees.md            the licence of each redistributed file,
                                     shipped as LICENSE-trees.md

Standard-library Python.

    python scripts/tree_licences.py          # print the classification
    python scripts/tree_licences.py --write  # and write the outputs (the
                                             # rebuild writes them too)
"""
import csv
import hashlib
import io
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))

LICENCES = ("Data", "tree_licences.csv")
WITHHELD = ("Data", "withheld_trees.csv")
TIPS_DIR = ("Data", "withheld_tips")
STATEMENT = ("Data", "LICENSE-trees.md")
RELEASE_TIPS = "registers/withheld_tips/"      # + release path + ".txt.gz"

REDISTRIBUTED = ("open", "open-share-alike")
WITHHELD_CLASSES = ("no-licence-stated",)

WITHHELD_COLUMNS = (
    "release_path", "paper_path", "kind", "study", "source_url", "licence",
    "redistribution", "licence_evidence_url", "sha256", "tips", "distinct_tips",
    "root_height", "roottotip_cv", "tips_list", "derivation")

STD = "standardized/trees/"
SRC = "standardized/canonical_sources/"
SUP = "standardized/superseded_canonicals/"
SITE = "site/data/partitions/"
SUP_DIR = ("Data", "_archive", "2026-10-revision", "superseded_canonicals", "trees")


def P(*parts):
    return os.path.join(ROOT, *parts)


def read_csv(*parts):
    with open(P(*parts), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _norm_ws(s):
    return re.sub(r"\s+", " ", s or "").strip()


# ---------------------------------------------------------------- trees ---
def leaf_depths(text):
    """Leaves of a Newick string as (label, root-to-tip depth) pairs.  The
    same reading as scripts/reproduce_paper.py of the release: whitespace
    outside quoted labels is ignored, the root's own branch length counts, a
    NaN branch length contributes zero depth."""
    s, quote = [], None
    for ch in text.strip():
        if quote:
            s.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            s.append(ch)
        elif not ch.isspace():
            s.append(ch)
    s = "".join(s)
    if s.endswith(";"):
        s = s[:-1]
    parent, bl, name, nkids = [-1], [0.0], [""], [0]
    cur, i, n = 0, 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            parent.append(cur); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[cur] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ",":
            p = parent[cur]
            parent.append(p); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[p] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ")":
            cur = parent[cur]
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),":
                j += 1
            lab, _, b = s[i:j].partition(":")
            name[cur] = lab.strip()
            b = b.strip()
            try:
                v = float(b) if b else 0.0
            except ValueError:
                v = 0.0
            if v != v:
                v = 0.0
            bl[cur] = v
            i = j
    depth = [0.0] * len(parent)
    for k in range(len(parent)):
        depth[k] = (depth[parent[k]] if parent[k] >= 0 else 0.0) + bl[k]
    return [(name[k], depth[k]) for k in range(len(parent)) if nkids[k] == 0]


def cv(xs):
    if not xs:
        return None
    m = sum(xs) / len(xs)
    if m == 0:
        return None
    return math.sqrt(sum((x - m) ** 2 for x in xs) / len(xs)) / abs(m)


def tip_sort_key(label):
    return (0, int(label), "") if label.lstrip("-").isdigit() else (1, 0, label)


def file_facts(path):
    data = open(path, "rb").read()
    leaves = leaf_depths(data.decode("utf-8", errors="replace"))
    depths = [d for _, d in leaves]
    c = cv(depths)
    return {
        "sha256": hashlib.sha256(data).hexdigest(),
        "tips": len(leaves),
        "distinct_tips": len({l for l, _ in leaves}),
        "root_height": "%.12g" % max(depths) if depths else "",
        "roottotip_cv": "%.6f" % c if c is not None else "",
        "labels": sorted((l for l, _ in leaves), key=tip_sort_key),
    }


# ----------------------------------------------------------- candidates ---
def candidates():
    """[(release path, paper path, kind)] of every tree file the release
    would carry if every licence allowed it."""
    out = []
    for r in read_csv("Data", "standardized", "metadata.csv"):
        out.append((STD + r["filename"], "Data/standardized/trees/" + r["filename"],
                    "catalogue tree"))
    for r in read_csv("Data", "canonical_sources.csv"):
        out.append((SRC + os.path.basename(r["source"]), r["source"],
                    "unpruned source of a reference tree"))
    for f in sorted(os.listdir(P(*SUP_DIR))):
        if f.endswith(".nwk"):
            out.append((SUP + f, "/".join(SUP_DIR + (f,)), "former reference tree"))
    base = P("Data", "site", "data", "partitions")
    for d, dirs, files in os.walk(base):
        dirs.sort()
        for f in sorted(files):
            if not f.endswith(".nwk"):
                continue
            rel_ = os.path.relpath(os.path.join(d, f), P("Data")).replace(os.sep, "/")
            if "/sub-phylos/" in rel_:
                out.append((rel_, "Data/" + rel_, "sub-clade tree of the browser site"))
    seen = set()
    for rp, _, _ in out:
        if rp in seen:
            raise SystemExit("tree_licences: %s is a candidate twice" % rp)
        seen.add(rp)
    return out


def licence_rows():
    rows = read_csv(*LICENCES)
    by_file = {}
    for r in rows:
        for f in (x.strip() for x in r["filenames"].split(";")):
            if not f or f.endswith("*"):
                continue
            if f in by_file:
                raise SystemExit("tree_licences.csv lists %s twice" % f)
            by_file[f] = r
    return rows, by_file


def classify():
    """[(release path, paper path, kind, licence row, how the row was found)]."""
    _, by_file = licence_rows()
    src_of = {os.path.basename(r["source"]): r["filename"]
              for r in read_csv("Data", "canonical_sources.csv")}
    out = []
    for rp, wp, kind in candidates():
        base = os.path.basename(rp)
        if rp.startswith(STD):
            key, how = base, "row of %s" % base
        elif rp.startswith(SRC):
            key = src_of[base]
            how = "row of %s, the reference tree built from it" % key
        elif rp.startswith(SUP):
            key, how = rp, "own row"
        else:
            if rp in by_file:
                key, how = rp, "own row"
            else:
                slug = rp.split("/sub-phylos/")[1].split("/")[0]
                std = P("Data", "standardized", "trees", slug + ".nwk")
                if not (os.path.exists(std) and open(P(*wp.split("/")), "rb").read()
                        == open(std, "rb").read()):
                    raise SystemExit("tree_licences: %s has no row in "
                                     "Data/tree_licences.csv and is not a copy of a "
                                     "catalogue tree" % rp)
                key, how = slug + ".nwk", "copy of %s%s.nwk" % (STD, slug)
        row = by_file.get(key)
        if row is None:
            raise SystemExit("tree_licences: no row of Data/tree_licences.csv covers "
                             "%s (looked for %s)" % (rp, key))
        cls = row["redistribution"]
        if cls not in REDISTRIBUTED + WITHHELD_CLASSES:
            raise SystemExit("tree_licences: %s has redistribution class %r"
                             % (rp, cls))
        out.append((rp, wp, kind, row, how))
    return out


# ------------------------------------------------------------ derivation ---
def _study_of():
    meta = {r["filename"]: r["study"] for r in read_csv("Data", "standardized", "metadata.csv")}
    return meta


def derivation(rp, wp, how):
    """How the file the paper used was derived from the source, from the
    registers that record it."""
    refs = {r["filename"]: r for r in read_csv("Data", "standardized", "tree_references.csv")}
    srcs = {os.path.basename(r["source"]): r for r in read_csv("Data", "canonical_sources.csv")}
    prep = {r["filename"]: r for r in read_csv("Data", "canonical_source_prep.csv")}
    base = os.path.basename(rp)
    parts = []
    if rp.startswith(STD):
        r = refs.get(base)
        if r:
            parts.append("Data/standardized/tree_references.csv: rule %s; %s"
                         % (r["rule"], _norm_ws(r["evidence"])))
        for s in srcs.values():
            if s["filename"] == base:
                parts.append("built from %s%s, whose origin is %s"
                             % (SRC, os.path.basename(s["source"]), s["origin"]))
                pr = prep.get(base)
                if pr and s["origin"] == pr["raw_source"]:
                    parts.append("Data/canonical_source_prep.csv: tree %s of %s, "
                                 "branch lengths x %s" % (pr["tree_index"],
                                                         pr["raw_source"],
                                                         pr["branch_scale"]))
    elif rp.startswith(SRC):
        r = srcs[base]
        parts.append("Data/canonical_sources.csv: source of %s, origin %s; %s"
                     % (r["filename"], r["origin"], _norm_ws(r["how"])))
        pr = prep.get(r["filename"])
        if pr and r["origin"] == pr["raw_source"]:
            parts.append("Data/canonical_source_prep.csv: tree %s of %s, branch "
                         "lengths x %s; %s" % (pr["tree_index"], pr["raw_source"],
                                               pr["branch_scale"], _norm_ws(pr["evidence"])))
        r2 = refs.get(wp)
        if r2:
            parts.append("Data/standardized/tree_references.csv: rule %s against %s"
                         % (r2["rule"], r2["reference"]))
    elif rp.startswith(SUP):
        r = refs.get(wp)
        if r:
            parts.append("Data/standardized/tree_references.csv: rule %s against %s; %s"
                         % (r["rule"], r["reference"], _norm_ws(r["evidence"])))
    else:
        slug = rp.split("/sub-phylos/")[1].split("/")[0]
        std = P("Data", "standardized", "trees", slug + ".nwk")
        if os.path.exists(std) and open(std, "rb").read() == \
                open(P(*wp.split("/")), "rb").read():
            parts.append("byte-identical copy of %s%s.nwk" % (STD, slug))
        else:
            info = os.path.join(os.path.dirname(P(*wp.split("/"))), "info.yaml")
            text = open(info, encoding="utf-8").read() if os.path.exists(info) else ""
            m = re.search(r"derivation: \|\n((?:    .*\n)+)", text)
            if m:
                parts.append(_norm_ws(m.group(1)))
            m = re.search(r"    note: \"([^\"]*)\"", text)
            if m:
                parts.append(_norm_ws(m.group(1)))
    if not parts:
        raise SystemExit("tree_licences: no register records how %s was derived" % rp)
    return " | ".join(parts)


# --------------------------------------------------------------- outputs ---
def tips_path(rp):
    return "/".join(TIPS_DIR) + "/" + rp + ".txt"


def build():
    """{'withheld': [row dicts], 'redistributed': [(rp, row)], 'tips': {paper
    path of the tips file: text}, 'files': {paper path: text}}."""
    study = _study_of()
    withheld, redistributed, tips = [], [], {}
    for rp, wp, kind, row, how in classify():
        if row["redistribution"] in REDISTRIBUTED:
            redistributed.append((rp, row))
            continue
        f = file_facts(P(*wp.split("/")))
        base = os.path.basename(rp)
        st = study.get(base) if not rp.startswith(SITE) else None
        if rp.startswith(SRC):
            st = study.get({os.path.basename(r["source"]): r["filename"] for r in
                            read_csv("Data", "canonical_sources.csv")}[base])
        withheld.append({
            "release_path": rp, "paper_path": wp, "kind": kind,
            "study": st or row["source_key"],
            "source_url": row["archive_url"] or row["licence_evidence_url"],
            "licence": row["licence"], "redistribution": row["redistribution"],
            "licence_evidence_url": row["licence_evidence_url"],
            "sha256": f["sha256"], "tips": str(f["tips"]),
            "distinct_tips": str(f["distinct_tips"]), "root_height": f["root_height"],
            "roottotip_cv": f["roottotip_cv"],
            "tips_list": RELEASE_TIPS + rp + ".txt.gz",
            "derivation": derivation(rp, wp, how),
        })
        if not row["archive_url"] and not row["licence_evidence_url"]:
            raise SystemExit("tree_licences: withheld %s has no source URL" % rp)
        tips[tips_path(rp)] = "\n".join(f["labels"]) + "\n"
    out = io.StringIO()
    w = csv.DictWriter(out, fieldnames=WITHHELD_COLUMNS, lineterminator="\n")
    w.writeheader()
    w.writerows(withheld)
    files = {"/".join(WITHHELD): out.getvalue(),
             "/".join(STATEMENT): statement(redistributed, withheld)}
    return {"withheld": withheld, "redistributed": redistributed, "tips": tips,
            "files": files}


def counts(b):
    """The counts the paper and the release print."""
    w = [r for r in b["withheld"] if r["release_path"].startswith(STD)]
    n_std = len(read_csv("Data", "standardized", "metadata.csv"))
    return {
        "n_catalogue": n_std,
        "n_redistributed": n_std - len(w),
        "n_withheld": len(w),
        "n_withheld_files": len(b["withheld"]),
        "n_redistributed_files": len(b["redistributed"]),
        "withheld_trees": sorted(os.path.basename(r["release_path"]) for r in w),
    }


def _md(s):
    return (s or "").replace("|", "\\|").replace("\n", " ")


LICENCE_TEXT = {
    "CC0-1.0": "[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)",
    "CC-BY-4.0": "[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)",
    "CC-BY-2.0": "[CC BY 2.0](https://creativecommons.org/licenses/by/2.0/)",
    "CC-BY-SA-4.0": "[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)",
    "GPL-3.0": "[GPL-3.0](LICENSES/GPL-3.0.txt)",
    "public domain (TreeBASE policy)": "public domain (TreeBASE policy)",
}


def statement(redistributed, withheld):
    by_lic = {}
    for rp, row in redistributed:
        by_lic.setdefault(row["release_licence"], []).append(rp)
    L = []
    L.append("# Licences of the tree files")
    L.append("")
    L.append("Generated by `scripts/tree_licences.py` of the paper repository from "
             "`registers/tree_licences.csv`, which records for each source the licence "
             "its archive states, the evidence (URL and quotation) and the date it was "
             "checked. Do not edit by hand.")
    L.append("")
    L.append("The catalogue's own work (the registers, the name layer, the code, the "
             "paper's generated files and the browser site) is under CC BY 4.0 "
             "(`LICENSE`). A tree file is the work of the study that published it, "
             "relabelled with the catalogue's tip identifiers and, for a reference tree, "
             "pruned to its group with one tip per species; each is redistributed under "
             "the licence of its source, listed below, and only where that licence "
             "allows redistribution. The changes made to every file are recorded in "
             "`standardized/tree_references.csv`, `registers/canonical_sources.csv` and "
             "`registers/canonical_source_prep.csv`.")
    L.append("")
    L.append("| Licence | Files |")
    L.append("|---|---|")
    for lic in sorted(by_lic):
        L.append("| %s | %d |" % (LICENCE_TEXT.get(lic, lic), len(by_lic[lic])))
    L.append("| not redistributed (no licence stated by the source) | %d |" % len(withheld))
    L.append("")
    sa = [rp for rp, row in redistributed if row["release_licence"] == "CC-BY-SA-4.0"]
    gpl = [rp for rp, row in redistributed if row["release_licence"] == "GPL-3.0"]
    if sa:
        L.append("**Share-alike.** %s %s distributed under CC BY-SA 4.0, the licence of "
                 "the Genome Taxonomy Database: the files as adapted here (tips relabelled, "
                 "off-group tips removed) are licensed under CC BY-SA 4.0, not under the "
                 "catalogue's CC BY 4.0, and any further adaptation must carry the same "
                 "licence." % (", ".join("`%s`" % x for x in sa),
                               "is" if len(sa) == 1 else "are"))
        L.append("")
    if gpl:
        L.append("**GPL.** %s %s the tree file of the GitHub repository named below, "
                 "which is licensed under the GNU General Public License v3.0; the "
                 "licence text is in `LICENSES/GPL-3.0.txt`. The file as distributed here "
                 "is modified (tip labels replaced by the catalogue's identifiers, branch "
                 "lengths written without trailing zeros) and is licensed under GPL-3.0."
                 % (", ".join("`%s`" % x for x in gpl), "is" if len(gpl) == 1 else "are"))
        L.append("")
    if withheld:
        L.append("**Not redistributed.** The sources of %d files state no licence, so "
                 "they are not in this release. `registers/withheld_trees.csv` lists each "
                 "with the source to download the original from, the sha256 and tip "
                 "count of the file the paper used and how it was derived; "
                 "`registers/withheld_tips/` lists the tip labels of each (species, not "
                 "trees). `scripts/reproduce_paper.py` checks a withheld file in full when "
                 "it is placed at its path." % len(withheld))
        L.append("")
        L.append("| File | Source | Licence stated |")
        L.append("|---|---|---|")
        for r in withheld:
            L.append("| `%s` | %s | %s |" % (r["release_path"], _md(r["source_url"]),
                                            _md(r["licence"])))
        L.append("")
    L.append("## Each redistributed file")
    L.append("")
    L.append("| File | Licence | Source | Attribution and conditions | Evidence |")
    L.append("|---|---|---|---|---|")
    for rp, row in sorted(redistributed):
        L.append("| `%s` | %s | %s | %s | %s |" % (
            rp, LICENCE_TEXT.get(row["release_licence"], row["release_licence"]),
            _md(row["archive_url"] or row["source_key"]), _md(row["conditions"]),
            _md(row["licence_evidence_url"])))
    return "\n".join(L) + "\n"


def write(b):
    """Write the outputs (scripts/rebuild_atlas.py writes the same files); a
    tip list of a file no longer withheld is removed."""
    out = dict(b["files"])
    out.update(b["tips"])
    for k, v in out.items():
        p = P(*k.split("/"))
        os.makedirs(os.path.dirname(p), exist_ok=True)
        old = open(p, encoding="utf-8").read() if os.path.exists(p) else None
        if old != v:
            with open(p, "w", encoding="utf-8", newline="") as fh:
                fh.write(v)
            print("wrote %s" % k)
    for d, _, fs in os.walk(P(*TIPS_DIR)):
        for f in fs:
            k = os.path.relpath(os.path.join(d, f), ROOT).replace(os.sep, "/")
            if k not in out:
                os.remove(os.path.join(d, f))
                print("removed %s" % k)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    b = build()
    c = counts(b)
    for r in b["withheld"]:
        print("withheld  %-75s %s" % (r["release_path"], r["source_url"]))
    print(c)
    if "--write" in argv:
        write(b)
    return 0


if __name__ == "__main__":
    sys.exit(main())
