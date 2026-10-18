#!/usr/bin/env python3
"""watch.py -- the monthly watch of the living atlas.

Finds phylogenetic trees deposited in a date window, measures each one with
the atlas's own rules, and screens it against the reference tree of its group.
It changes nothing in the atlas: it writes a run folder with a ledger of every
tree seen and a report of the candidates, which a curation step then reads.

STEPS
  harvest   records of the window from Zenodo, Figshare and Dryad whose text
            matches a phylogeny query; the tree files of each record (by file
            name), downloaded (Dryad through its API account: an OAuth token
            from ~/.config/atlas-watch/dryad.json or the vault; without them
            Dryad records are listed, not measured).
  measure   the first tree of each file (Newick or Nexus, inside zip, tar or
            gzip archives too): tips, species keys by the atlas rule
            (scripts/canonical_trees.py species_key, D4/D13), number of trees
            in the file, root-to-tip coefficient of variation.
  assign    each genus is classified in the Catalogue of Life COL26.9 base
            release (ChecklistBank dataset 316321, the release of the paper)
            and each species key is assigned to the group whose definition in
            Data/col_denominators.csv (col_taxa: "A + B - C") contains it; the
            tree's group is the group of most of its species.
  screen    against the two reference trees of each group (decision D65), read
            from Data/coverage_basis.csv: the secondary tree, the largest tree
            of any kind (the reference tree), and the primary tree, the largest
            tree dated in absolute time (the reference tree where its
            time_scale is absolute).  A tree is a candidate when it is
            species-level (at least MIN_SPECIES_SHARE of its tips carry a
            species key), most of its species are in one group (at least
            MIN_GROUP_SHARE), and
              LARGER        its species in the group are at least SCREEN_RATIO
                            times the secondary tree's accepted names;
              FIRST         the group has no reference tree and the tree holds
                            at least MIN_FIRST_SPECIES species of it;
              DATED-LARGER  it is ultrametric (CV <= CHRONOGRAM_TOLERANCE) and
                            at least SCREEN_RATIO times the primary tree;
              DATED-FIRST   it is ultrametric, the group has no primary tree,
                            and it holds at least MIN_FIRST_SPECIES species;
              SET           an ultrametric tree set as large as the primary tree
                            where that archives no date uncertainty;
              REDEPOSIT?    added when at least 95% of the tree's species in
                            the group are in the reference tree and fewer than
                            5% are new: most likely the reference tree itself,
                            deposited again by a later study.
            Species keys are not yet resolved to accepted names, so the screen
            ratio is below one: curation counts accepted names by
            scripts/match_col.py before any change.

USAGE
  python living/watch.py run --since 2026-06-01 [--until 2026-10-14]
  python living/watch.py run --last-month          # the previous calendar month
Output: living/runs/<since>_<until>/ (ledger.csv.gz, records.csv.gz, report.md);
downloads and the genus cache in living/cache/ (not tracked).
Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import http.client as httpclient
import io
import json
import math
import os
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import Counter, OrderedDict
from concurrent.futures import ThreadPoolExecutor
import threading

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from canonical_trees import parse, leaves, species_key, _blf  # noqa: E402

DATA = os.path.join(ROOT, "Data")
CACHE = os.path.join(HERE, "cache")
RUNS = os.path.join(HERE, "runs")

COL_DATASET = 316321                 # COL26.9 base release, as in the paper
CHRONOGRAM_TOLERANCE = 0.05          # rebuild_atlas.py, \RuleChronCV
MIN_SPECIES_SHARE = 0.5
MIN_GROUP_SHARE = 0.6
SCREEN_RATIO = 0.9
MIN_FIRST_SPECIES = 20
MAX_DOWNLOAD = 300 * 1024 * 1024     # bytes per tree file
MAX_ARCHIVE = 50 * 1024 * 1024       # bytes per archive (zip, tar, gz)
MAX_TREE_TEXT = 400 * 1024 * 1024    # characters read from one tree file
QUERY_TERMS = ["phylogeny", "phylogenetic", "phylogenomic", "phylogenomics",
               "chronogram", "time-calibrated", "time-tree", "timetree",
               "divergence times"]
TREE_EXT = (".tre", ".tree", ".trees", ".nwk", ".newick", ".nex", ".nexus",
            ".treefile", ".contree", ".phy.treefile", ".tree.txt", ".nhx")
ARCHIVE_EXT = (".zip", ".tar.gz", ".tgz", ".tar", ".gz")
UA = "phylo-atlas-watch/1.0 (mailto:richtf@usi.ch)"


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------

def http(url, data=None, headers=None, tries=4, raw=False, limit=None,
         timeout=120):
    h = {"User-Agent": UA, "Accept": "application/json"}
    h.update(headers or {})
    body = json.dumps(data).encode() if data is not None else None
    if body is not None:
        h["Content-Type"] = "application/json"
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=body, headers=h)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if raw:
                    return r.read(limit) if limit else r.read()
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code in (400, 401, 403, 404, 410):
                raise
            time.sleep(2 ** i * 3)
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                httpclient.HTTPException, ValueError):
            time.sleep(2 ** i * 3)
    raise RuntimeError("giving up on %s" % url)


def is_tree_file(name):
    n = name.lower()
    if n.endswith(TREE_EXT):
        return True
    if n.endswith(".txt") and re.search(r"tree|phylo|chrono", n):
        return True
    return False


def is_archive(name):
    return name.lower().endswith(ARCHIVE_EXT)


def textual_hit(*texts):
    t = " ".join(x or "" for x in texts).lower()
    return any(w in t for w in ("phylogen", "chronogram", "time-calibrated",
                                "time tree", "timetree", "divergence time",
                                "species tree", "dated tree"))


# ---------------------------------------------------------------------------
# Harvest
# ---------------------------------------------------------------------------

def _zenodo_record(h):
    md = h.get("metadata", {})
    files = [{"name": f.get("key", ""), "size": f.get("size", 0),
              "url": f.get("links", {}).get("self", "")}
             for f in h.get("files", [])]
    return OrderedDict(
        source="zenodo", record=str(h.get("id")),
        doi=h.get("doi", ""), date=md.get("publication_date", ""),
        title=re.sub(r"\s+", " ", md.get("title", "")),
        related=";".join(r.get("identifier", "") for r in
                         md.get("related_identifiers", [])
                         if r.get("relation") in
                         ("isSupplementTo", "isCitedBy", "isPartOf")),
        files=files, downloadable=True,
        kind=(md.get("resource_type") or {}).get("type", ""),
        text_hit=textual_hit(md.get("title"), md.get("description")))


def harvest_zenodo(since, until):
    """In 7-day slices, so no query pages deep (deep pages fail), 25 records
    a page (the limit without a token).  Each slice's record count is checked
    against the total Zenodo reports; a shortfall is returned as an _error
    row, which makes the run incomplete."""
    terms = " OR ".join('"%s"' % t if " " in t or "-" in t else t
                        for t in QUERY_TERMS)
    out = []
    d0, d1 = dt.date.fromisoformat(since), dt.date.fromisoformat(until)
    while d0 <= d1:
        e = min(d0 + dt.timedelta(days=6), d1)
        q = "(%s) AND publication_date:[%s TO %s]" % (terms, d0, e)
        got, total, page = [], None, 1
        while True:
            url = ("https://zenodo.org/api/records?" + urllib.parse.urlencode(
                {"q": q, "size": 25, "page": page, "sort": "mostrecent",
                 "all_versions": "false"}))
            d = None
            for attempt in range(3):          # a page that fails is retried
                try:
                    d = http(url)
                    break
                except Exception as err:
                    last = err
                    time.sleep(60)
            if d is None:
                out.append({"_error": "zenodo %s..%s page %d: %s" %
                            (d0, e, page, last)})
                break
            total = d["hits"].get("total", total)
            hits = d["hits"]["hits"]
            got.extend(_zenodo_record(h) for h in hits)
            if len(hits) < 25 or page >= 400:
                break
            page += 1
        if total is not None and len(got) < total:
            out.append({"_error": "zenodo %s..%s: %d of %d records" %
                        (d0, e, len(got), total)})
        out.extend(got)
        d0 = e + dt.timedelta(days=1)
    return out


def harvest_figshare(since, until):
    out, seen = [], set()
    for term in ("phylogeny", "phylogenetic tree", "phylogenomic",
                 "time-calibrated", "chronogram"):
        page = 1
        while True:
            d = http("https://api.figshare.com/v2/articles/search",
                     data={"search_for": term, "published_since": since,
                           "page_size": 100, "page": page})
            for a in d:
                if a["id"] in seen or a.get("published_date", "")[:10] > until:
                    continue
                seen.add(a["id"])
                out.append(OrderedDict(
                    source="figshare", record=str(a["id"]),
                    doi=a.get("doi", ""),
                    date=a.get("published_date", "")[:10],
                    title=re.sub(r"\s+", " ", a.get("title", "")),
                    related="", files=None, downloadable=True,
                    text_hit=textual_hit(a.get("title"))))
            if len(d) < 100 or page >= 20:
                break
            page += 1
    return out


def figshare_files(rec):
    d = http("https://api.figshare.com/v2/articles/%s" % rec["record"])
    rec["related"] = d.get("resource_doi", "") or ""
    rec["text_hit"] = rec["text_hit"] or textual_hit(d.get("description"))
    return [{"name": f.get("name", ""), "size": f.get("size", 0),
             "url": f.get("download_url", "")} for f in d.get("files", [])]


def harvest_dryad(since, until):
    out, seen = [], set()
    for term in QUERY_TERMS:
        page = 1
        while True:
            url = ("https://datadryad.org/api/v2/search?" + urllib.parse.urlencode(
                {"q": term, "per_page": 100, "page": page,
                 "publishedSince": since}))
            d = http(url)
            ds = d.get("_embedded", {}).get("stash:datasets", [])
            for x in ds:
                ident = x.get("identifier", "")
                date = x.get("publicationDate", "") or ""
                if ident in seen or not date or date > until or date < since:
                    continue
                seen.add(ident)
                rel = ";".join(r.get("identifier", "") for r in
                               x.get("relatedWorks", []) or []
                               if r.get("relationship") == "primary_article")
                out.append(OrderedDict(
                    source="dryad", record=ident, doi=ident.replace("doi:", ""),
                    date=date, title=re.sub(r"<[^>]+>|\s+", " ",
                                            x.get("title", "")).strip(),
                    related=rel, files=None,
                    downloadable=Dryad.credentials() is not None,
                    version=x.get("_links", {}).get("stash:version", {})
                    .get("href", ""),
                    text_hit=textual_hit(x.get("title"), x.get("abstract"))))
            if len(ds) < 100 or page >= 20:
                break
            page += 1
    return out


class Dryad(object):
    """OAuth client-credentials token for Dryad file downloads.  Credentials
    from DRYAD_CLIENT_ID / DRYAD_CLIENT_SECRET, else
    ~/.config/atlas-watch/dryad.json, else `vault get dryad` (the Macs)."""
    _token, _when, _cred = None, 0.0, False

    @classmethod
    def credentials(cls):
        if cls._cred is not False:
            return cls._cred
        cid, sec = os.environ.get("DRYAD_CLIENT_ID"), os.environ.get("DRYAD_CLIENT_SECRET")
        f = os.path.expanduser("~/.config/atlas-watch/dryad.json")
        if cid and sec:
            cls._cred = (cid, sec)
        elif os.path.exists(f):
            d = json.load(open(f))
            cls._cred = (d["client_id"], d["client_secret"])
        else:
            try:
                import subprocess
                g = lambda k: subprocess.check_output(
                    ["vault", "get", "dryad", k], text=True,
                    stderr=subprocess.DEVNULL).strip()
                cls._cred = (g("username"), g("password"))
            except Exception:
                cls._cred = None
        return cls._cred

    @classmethod
    def token(cls):
        if cls._token and time.time() - cls._when < 9 * 3600:
            return cls._token
        cred = cls.credentials()
        body = urllib.parse.urlencode({"grant_type": "client_credentials",
                                       "client_id": cred[0],
                                       "client_secret": cred[1]}).encode()
        req = urllib.request.Request("https://datadryad.org/oauth/token", data=body,
                                     headers={"Accept": "application/json",
                                              "User-Agent": UA})
        cls._token = json.loads(urllib.request.urlopen(req, timeout=60).read())[
            "access_token"]
        cls._when = time.time()
        return cls._token


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


_DRYAD_LOCK = threading.Lock()


def dryad_download(url, timeout=60, tries=6):
    """One Dryad download at a time (parallel requests draw 429), waiting what
    Retry-After asks on a 429."""
    with _DRYAD_LOCK:
        for i in range(tries):
            try:
                return _dryad_get(url, timeout)
            except urllib.error.HTTPError as e:
                if e.code != 429 or i == tries - 1:
                    raise
                wait = e.headers.get("Retry-After") or ""
                time.sleep(int(wait) if wait.isdigit() else 30 * (i + 1))


def _dryad_get(url, timeout):
    """Dryad answers 302 to a signed S3 URL, which rejects the bearer token:
    the redirect is followed without it."""
    req = urllib.request.Request(url, headers={
        "Authorization": "Bearer " + Dryad.token(), "User-Agent": UA})
    try:
        with urllib.request.build_opener(_NoRedirect).open(req, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError as e:
        if e.code not in (301, 302, 303, 307, 308):
            raise
        with urllib.request.urlopen(e.headers["Location"], timeout=timeout) as r:
            return r.read()


def dryad_files(rec):
    if not rec.get("version"):
        return []
    d = http("https://datadryad.org%s/files" % rec["version"])
    return [{"name": f.get("path", ""), "size": f.get("size", 0),
             "url": "https://datadryad.org" + f.get("_links", {})
             .get("stash:download", {}).get("href", "")}
            for f in d.get("_embedded", {}).get("stash:files", [])
            if f.get("_links", {}).get("stash:download")]


# ---------------------------------------------------------------------------
# Reading trees
# ---------------------------------------------------------------------------

NEXUS_TREE_RE = re.compile(r"^\s*u?tree\s+[^=]+=\s*(?:\[[^\]]*\]\s*)*", re.I | re.M)
# node intervals as BEAST/TreeAnnotator (height_95%_HPD) or MCMCTree ([&95%={..}]) write them
INTERVAL_RE = re.compile(r"95%_?HPD|\[&95%=|_95%=", re.I)


def tree_texts(raw):
    """(first tree as Newick with translated labels, number of trees)."""
    text = raw.decode("utf-8", errors="replace")
    if re.match(r"\s*#nexus", text, re.I):
        trans = {}
        m = re.search(r"translate\s+(.*?);", text, re.I | re.S)
        if m:
            for part in re.split(r",\s*(?=\S)", m.group(1)):
                bits = part.strip().split(None, 1)
                if len(bits) == 2:
                    trans[bits[0]] = bits[1].strip().strip(",")
        starts = [x.end() for x in NEXUS_TREE_RE.finditer(text)]
        if not starts:
            return None, 0
        first = text[starts[0]:text.index(";", starts[0]) + 1]
        return (first, trans), len(starts)
    # Newick: statements that start with '('
    stmts = [s for s in re.split(r";\s*", text) if s.strip().startswith("(")]
    if not stmts:
        return None, 0
    return (stmts[0] + ";", {}), len(stmts)


def measure_tree(raw):
    got, n_trees = tree_texts(raw)
    if not got:
        return None
    nwk, trans = got
    try:
        root = parse(nwk)
    except Exception:
        return None
    lv = leaves(root)
    if len(lv) < 4:
        return None
    labels = []
    for leaf in lv:
        lab = leaf.label
        if trans and lab in trans:
            lab = trans[lab].strip("'\"")
        labels.append(lab)
    keys = set()
    n_key = 0
    for i, lab in enumerate(labels):
        k, ok = species_key(lab.replace("'", ""), str(i))
        if ok:
            keys.add(k)
            n_key += 1
    # root-to-tip coefficient of variation (only when every branch has a length)
    cv = ""
    has_bl = all(n.bl is not None for n in lv)
    if has_bl:
        depth = {}
        stack = [(root, 0.0)]
        ok = True
        tips = []
        while stack:
            n, d = stack.pop()
            dd = d + (_blf(n.bl) if n is not root else 0.0)
            if math.isnan(dd):
                ok = False
                break
            if n.kids:
                stack.extend((k, dd) for k in n.kids)
            else:
                tips.append(dd)
        if ok and tips:
            mean = sum(tips) / len(tips)
            if mean > 0:
                var = sum((x - mean) ** 2 for x in tips) / len(tips)
                cv = math.sqrt(var) / mean
    return {"tips": len(lv), "keyed_tips": n_key, "species": sorted(keys),
            "n_trees": n_trees, "cv": cv,
            "intervals": "yes" if INTERVAL_RE.search(nwk) else "",
            "root_height": "" if cv == "" else max(tips)}


def files_in(name, blob):
    """(member name, bytes) of every tree file in a download."""
    n = name.lower()
    try:
        if n.endswith(".zip"):
            z = zipfile.ZipFile(io.BytesIO(blob))
            for info in z.infolist():
                if is_tree_file(info.filename) and info.file_size <= MAX_TREE_TEXT:
                    yield name + "!" + info.filename, z.read(info)
        elif n.endswith((".tar.gz", ".tgz", ".tar")):
            t = tarfile.open(fileobj=io.BytesIO(blob))
            for m in t.getmembers():
                if m.isfile() and is_tree_file(m.name) and m.size <= MAX_TREE_TEXT:
                    yield name + "!" + m.name, t.extractfile(m).read()
        elif n.endswith(".gz"):
            inner = name[:-3]
            if is_tree_file(inner):
                yield name, gzip.decompress(blob)[:MAX_TREE_TEXT]
        elif is_tree_file(name):
            yield name, blob[:MAX_TREE_TEXT]
    except (zipfile.BadZipFile, tarfile.TarError, OSError, EOFError):
        return


# ---------------------------------------------------------------------------
# Groups
# ---------------------------------------------------------------------------

def load_groups():
    groups = []
    with open(os.path.join(DATA, "col_denominators.csv"), newline="") as fh:
        for r in csv.reader(fh):
            if not r or r[0].startswith("#") or r[0] == "partition":
                continue
            expr = r[1]
            plus, minus = set(), set()
            for sign, name in re.findall(r"(^|[+-])\s*([^+-]+?)\s*(?=[+-]|$)",
                                         expr):
                (minus if sign == "-" else plus).add(name.strip())
            groups.append((r[0], plus, minus))
    return groups


def load_basis():
    out = {}
    with open(os.path.join(DATA, "coverage_basis.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            out[r["partition"]] = r
    return out


def load_primary():
    """Groups whose primary tree (D65) is not their reference tree, from
    living/primary_trees.csv: {group: (tree, accepted species)}."""
    out = {}
    p = os.path.join(HERE, "primary_trees.csv")
    if os.path.exists(p):
        with open(p, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                out[r["group"]] = (r["tree"], int(r["accepted_species"]))
    return out


PRIMARY = load_primary()

_REF_KEYS = None


def reference_keys():
    """Species keys of the tips of each reference tree, by group
    (col_resolution.csv), and of every other shipped tree, by file name
    (col_resolution_noncanonical.csv, for primary trees that are not the
    reference tree)."""
    global _REF_KEYS
    if _REF_KEYS is None:
        _REF_KEYS = {}
        for fn, col in (("col_resolution.csv", "partition"),
                        ("col_resolution_noncanonical.csv", "filename")):
            path = os.path.join(DATA, "standardized", fn)
            if not os.path.exists(path):
                continue
            with open(path, newline="", encoding="utf-8") as fh:
                for r in csv.DictReader(fh):
                    if r["query_source"] == "species_key":
                        _REF_KEYS.setdefault(r[col], set()).add(r["query"])
    return _REF_KEYS


class GenusClassifier(object):
    """Genus -> set of names in its COL26.9 classification, cached."""

    def __init__(self):
        os.makedirs(CACHE, exist_ok=True)
        self.path = os.path.join(CACHE, "col269_genus_classification.json")
        self.cache = {}
        if os.path.exists(self.path):
            self.cache = {k: [n for n in v if n] for k, v in
                          json.load(open(self.path)).items()}
        self.dirty = 0

    def _match(self, q, rank):
        url = ("https://api.checklistbank.org/dataset/%d/match/nameusage?" %
               COL_DATASET + urllib.parse.urlencode({"q": q, "rank": rank}))
        try:
            u = http(url).get("usage") or {}
        except Exception:
            return []
        if not u:
            return []
        names = [u.get("name", "")] + [c.get("name", "") for c in
                                       u.get("classification", [])]
        return [n for n in names if n]

    def get(self, genus, examples=()):
        """A genus name shared by homonyms in COL does not match; it is then
        classified through up to three of the tree's species of that genus."""
        if genus in self.cache:             # an empty answer is final too
            return self.cache[genus]
        names = self._match(genus, "genus")
        for sp in list(examples)[:3]:
            if names:
                break
            names = self._match(sp.replace("_", " "), "species")
        self.cache[genus] = names
        self.dirty += 1
        return names

    def save(self):
        with open(self.path, "w") as fh:
            json.dump(self.cache, fh)


def group_of(names, groups):
    """names: the genus first, then its classification towards the root.
    Where two definitions contain the genus, the group whose taxon is
    nearest the genus wins."""
    s = set(names)
    hits = [(min(names.index(t) for t in plus if t in s), g)
            for g, plus, minus in groups if s & plus and not s & minus]
    return min(hits)[1] if hits else ""


# ---------------------------------------------------------------------------
# Screen
# ---------------------------------------------------------------------------

def screen(m, groups, basis, clf, max_genera=400):
    sp = m["species"]
    if not sp or m["keyed_tips"] < MIN_SPECIES_SHARE * m["tips"]:
        return dict(group="", group_species=0, group_share="", verdict="",
                    reason="not species-level")
    by_genus = Counter(k.split("_")[0] for k in sp)
    examples = {}
    for k in sp:
        examples.setdefault(k.split("_")[0], []).append(k)
    gen_group = {}
    for g, _ in by_genus.most_common(max_genera):
        gen_group[g] = group_of(clf.get(g, examples.get(g, ())), groups)
    counts, keys_of = Counter(), {}
    for k in sp:
        gg = gen_group.get(k.split("_")[0], "")
        if gg:
            counts[gg] += 1
            keys_of.setdefault(gg, set()).add(k)
    if not counts:
        return dict(group="", group_species=0, group_share="", verdict="",
                    reason="no genus classified")
    grp, n = counts.most_common(1)[0]
    share = n / float(sum(counts.values()))
    out = dict(group=grp, group_species=n, group_share=round(share, 3))
    if share < MIN_GROUP_SHARE:
        out.update(verdict="", reason="spans groups (%s)" %
                   ", ".join("%s %d" % kv for kv in counts.most_common(3)))
        return out
    b = basis.get(grp, {})
    ref = int(float(b.get("c0_units") or 0))
    ultra = m["cv"] != "" and m["cv"] <= CHRONOGRAM_TOLERANCE
    # D65: the secondary reference tree is the largest tree of any kind (the
    # reference tree of the release); the primary is the largest tree dated in
    # absolute time, which on the release is the reference tree when that is
    # an absolute chronogram.  An ultrametric tree is screened against both;
    # curation checks that its time scale is absolute.
    primary = b.get("canonical_tree") if b.get("time_scale") == "absolute" else ""
    pref = ref if primary else 0
    if grp in PRIMARY:                       # living/primary_trees.csv
        primary, pref = PRIMARY[grp]
    verdicts = []
    if not b.get("canonical_tree"):
        if n >= MIN_FIRST_SPECIES:
            verdicts.append("FIRST")
    elif n >= SCREEN_RATIO * ref:
        verdicts.append("LARGER")
    if ultra:
        if not primary:
            if n >= MIN_FIRST_SPECIES:
                verdicts.append("DATED-FIRST")
        elif n >= SCREEN_RATIO * pref:
            verdicts.append("DATED-LARGER")
        if m["n_trees"] > 1 and b.get("uncertainty_available") != "yes" and \
                n >= SCREEN_RATIO * pref:
            verdicts.append("SET")
    dated_only = verdicts and all(v.startswith(("DATED", "SET")) for v in verdicts)
    rk = reference_keys().get(primary if dated_only and grp in PRIMARY else grp, set())
    mine = keys_of[grp]
    overlap = len(mine & rk) / float(len(mine)) if rk else 0.0
    new = len(mine - rk)
    if verdicts and rk and overlap >= 0.95 and new < 0.05 * len(mine):
        verdicts.append("REDEPOSIT?")
    out.update(reference=b.get("canonical_tree", ""), reference_species=ref,
               primary_reference=primary, primary_species=pref,
               reference_overlap=round(overlap, 3) if rk else "",
               not_in_reference=new if rk else "",
               verdict="+".join(verdicts),
               reason="" if verdicts else "smaller than the reference tree"
               if b.get("canonical_tree") else "fewer than %d species" %
               MIN_FIRST_SPECIES)
    return out


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def compact(ledger, keep=20):
    """Trees that are not species-level are summarised per downloaded file
    when there are more than `keep` of them (gene-tree archives hold
    thousands): one row with the count replaces them."""
    by_file = Counter((r.get("source"), r.get("record"),
                       str(r.get("file", "")).split("!")[0])
                      for r in ledger if r.get("reason") == "not species-level")
    out, done = [], set()
    for r in ledger:
        k = (r.get("source"), r.get("record"), str(r.get("file", "")).split("!")[0])
        if r.get("reason") == "not species-level" and by_file[k] > keep:
            if k not in done:
                done.add(k)
                row = dict((c, r.get(c, "")) for c in
                           ("source", "record", "doi", "date", "title",
                            "related"))
                row.update(file=k[2], n_trees=by_file[k],
                           reason="%d trees not species-level" % by_file[k])
                out.append(row)
            continue
        out.append(r)
    return out


HARVEST_ERRORS = []

LEDGER_COLS = ["source", "record", "doi", "date", "title", "related", "file",
               "tips", "keyed_tips", "species", "n_trees", "cv", "ultrametric",
               "intervals",
               "group", "group_species", "group_share", "reference",
               "reference_species", "primary_reference", "primary_species",
               "reference_overlap", "not_in_reference",
               "verdict", "reason"]


# ---------------------------------------------------------------------------
# Seen records and the search of the whole archive
# ---------------------------------------------------------------------------

STATE = os.path.join(HERE, "state")
SEEN = os.path.join(STATE, "seen.csv.gz")
SEARCH_MODE = False


def load_seen():
    out = set()
    if os.path.exists(SEEN):
        with gzip.open(SEEN, "rt", newline="") as fh:
            for r in csv.DictReader(fh):
                out.add((r["source"], r["record"]))
    return out


def remember_seen(records, run_name):
    os.makedirs(STATE, exist_ok=True)
    rows = []
    if os.path.exists(SEEN):
        with gzip.open(SEEN, "rt", newline="") as fh:
            rows = list(csv.DictReader(fh))
    have = {(r["source"], r["record"]) for r in rows}
    for r in records:
        k = (r["source"], r["record"])
        if k not in have:
            have.add(k)
            rows.append({"source": k[0], "record": k[1], "first_seen": run_name})
    with gzip.open(SEEN, "wt", newline="") as fh:
        w = csv.DictWriter(fh, ["source", "record", "first_seen"])
        w.writeheader()
        w.writerows(rows)


def group_terms(groups, basis):
    """Each eukaryotic group's search terms: its Catalogue of Life taxa
    (col_denominators.csv, the '+' names) and its own name."""
    out = OrderedDict()
    for g, plus, minus in groups:
        b = basis.get(g, {})
        if b.get("is_prokaryotic") == "yes":
            continue
        names = sorted(plus)
        if not g.lower().startswith(("other", "minor")):
            names.append(g)
        out[g] = names
    return out


def _zenodo_query(q, pages):
    out = []
    for page in range(1, pages + 1):
        url = ("https://zenodo.org/api/records?" + urllib.parse.urlencode(
            {"q": q, "size": 25, "page": page, "sort": "bestmatch",
             "all_versions": "false"}))
        try:
            hits = http(url)["hits"]["hits"]
        except Exception as e:
            out.append({"_error": "zenodo search page %d: %s" % (page, e)})
            break
        out.extend(_zenodo_record(h) for h in hits)
        if len(hits) < 25:
            break
    return out


def search_zenodo(terms, pages=4):
    phylo = " OR ".join('"%s"' % t if " " in t or "-" in t else t
                        for t in QUERY_TERMS)
    taxa = " OR ".join('"%s"' % t if " " in t else t for t in terms)
    return _zenodo_query("(%s) AND (%s)" % (phylo, taxa), pages)


def search_dryad(terms, per_page=100):
    out, seen = [], set()
    for t in terms:
        url = ("https://datadryad.org/api/v2/search?" + urllib.parse.urlencode(
            {"q": "%s phylogeny" % t, "per_page": per_page, "page": 1}))
        try:
            ds = http(url).get("_embedded", {}).get("stash:datasets", [])
        except Exception as e:
            out.append({"_error": "dryad search %s: %s" % (t, e)})
            continue
        for x in ds:
            ident = x.get("identifier", "")
            if ident in seen:
                continue
            seen.add(ident)
            rel = ";".join(r.get("identifier", "") for r in
                           x.get("relatedWorks", []) or []
                           if r.get("relationship") == "primary_article")
            out.append(OrderedDict(
                source="dryad", record=ident, doi=ident.replace("doi:", ""),
                date=x.get("publicationDate", "") or "",
                title=re.sub(r"<[^>]+>|\s+", " ", x.get("title", "")).strip(),
                related=rel, files=None,
                downloadable=Dryad.credentials() is not None,
                version=x.get("_links", {}).get("stash:version", {}).get("href", ""),
                text_hit=textual_hit(x.get("title"), x.get("abstract"))))
    return out


def search_figshare(terms, page_size=100):
    out, seen = [], set()
    for t in terms:
        try:
            d = http("https://api.figshare.com/v2/articles/search",
                     data={"search_for": "%s phylogeny" % t, "page_size": page_size,
                           "page": 1})
        except Exception as e:
            out.append({"_error": "figshare search %s: %s" % (t, e)})
            continue
        for a in d:
            if a["id"] in seen:
                continue
            seen.add(a["id"])
            out.append(OrderedDict(
                source="figshare", record=str(a["id"]), doi=a.get("doi", ""),
                date=a.get("published_date", "")[:10],
                title=re.sub(r"\s+", " ", a.get("title", "")),
                related="", files=None, downloadable=True,
                text_hit=textual_hit(a.get("title"))))
    return out


def search_harvest(selected):
    """The harvest list of a search: every source, every selected group."""
    groups, basis = load_groups(), load_basis()
    terms = group_terms(groups, basis)
    sel = [g for g in terms if not selected or g in selected]

    def each(fn):
        def go():
            out, ids = [], set()
            for g in sel:
                for r in fn(terms[g]):
                    k = (r.get("source"), r.get("record"))
                    if "_error" in r or k not in ids:
                        ids.add(k)
                        out.append(r)
            return out
        return go
    return [("zenodo", each(search_zenodo)), ("figshare", each(search_figshare)),
            ("dryad", each(search_dryad))]


def run(since, until, out_dir, log, harvest=None):
    """harvest: [(source name, function returning records)]; by default the
    three date-window harvesters."""
    groups, basis, clf = load_groups(), load_basis(), GenusClassifier()
    dl_dir = os.path.join(CACHE, "downloads")
    os.makedirs(dl_dir, exist_ok=True)
    records = []
    if harvest is None:
        harvest = [(n, (lambda f: lambda: f(since, until))(f)) for n, f in
                   (("zenodo", harvest_zenodo), ("figshare", harvest_figshare),
                    ("dryad", harvest_dryad))]
    for name, fn in harvest:
        try:
            got = fn()
        except Exception as e:                       # one source down
            log("harvest %s failed: %s" % (name, e))
            HARVEST_ERRORS.append("%s: harvest failed: %s" % (name, e))
            got = []
        for e in [r for r in got if "_error" in r]:
            log("harvest %s incomplete: %s" % (name, e["_error"]))
            HARVEST_ERRORS.append("%s: %s" % (name, e["_error"]))
        got = [r for r in got if "_error" not in r]
        log("harvest %s: %d records" % (name, len(got)))
        records.extend(got)
    seen = load_seen()
    if harvest is not None and SEARCH_MODE:
        fresh = [r for r in records if (r["source"], r["record"]) not in seen]
        log("not seen before: %d of %d records" % (len(fresh), len(records)))
        records = fresh
    remember_seen(records, os.path.basename(out_dir))
    ledger, rec_rows, jobs = [], [], []
    for i, rec in enumerate(records):
        if i % 250 == 0:
            log("records %d/%d, trees measured %d" % (
                i, len(records), sum(1 for r in ledger if r.get("tips"))))
            clf.save()
        if rec["files"] is None and not rec["text_hit"]:
            rec["files"] = []        # title does not name a tree: not listed
        try:
            if rec["files"] is None:
                rec["files"] = (figshare_files(rec) if rec["source"] == "figshare"
                                else dryad_files(rec))
        except Exception as e:
            log("files of %s %s: %s" % (rec["source"], rec["record"], e))
            rec["files"] = []
        cand = [f for f in rec["files"]
                if is_tree_file(f["name"]) or (
                    is_archive(f["name"]) and rec["text_hit"] and
                    rec.get("kind") != "software" and not
                    re.search(r"matri|align|fasta|fastq|sequence|reads|raw|"
                              r"image|photo|video|figure", f["name"], re.I))]
        rec_rows.append(OrderedDict(
            (k, rec.get(k, "")) for k in ("source", "record", "doi", "date",
                                          "title", "related", "text_hit")))
        rec_rows[-1]["tree_files"] = ";".join(f["name"] for f in cand)
        rec_rows[-1]["downloadable"] = rec["downloadable"]
        if not cand:
            continue
        base = OrderedDict((k, rec.get(k, "")) for k in
                           ("source", "record", "doi", "date", "title", "related"))
        if not rec["downloadable"]:
            for f in cand:
                row = dict(base, file=f["name"], reason="not downloadable "
                           "(Dryad API account needed)")
                ledger.append(row)
            continue
        for f in cand:
            cap = MAX_ARCHIVE if is_archive(f["name"]) else MAX_DOWNLOAD
            if not f.get("url") or (f.get("size") or 0) > cap:
                ledger.append(dict(base, file=f["name"],
                                   reason="too large or no link (%s bytes)" %
                                   f.get("size")))
                continue
            local = os.path.join(dl_dir, "%s_%s_%s" % (
                rec["source"], re.sub(r"[^A-Za-z0-9._-]", "_", rec["record"]),
                re.sub(r"[^A-Za-z0-9._-]", "_", f["name"])[-120:]))
            jobs.append((base, f, local))
    log("downloading %d files" % len(jobs))

    def fetch(job):
        base, f, local = job
        if os.path.exists(local):
            return None
        try:
            if f["url"].startswith("https://datadryad.org"):
                blob = dryad_download(f["url"])
            else:
                blob = http(f["url"], raw=True, headers={"Accept": "*/*"},
                            tries=2, timeout=60)
        except Exception as e:
            return "download failed: %s" % e
        try:
            with open(local + ".part", "wb") as fh:
                fh.write(blob)
            os.replace(local + ".part", local)
        except OSError as e:
            return "could not store the download: %s" % e
        return None

    with ThreadPoolExecutor(6) as pool:
        errors = list(pool.map(fetch, jobs))
    log("measuring")
    measured = []
    for (base, f, local), err in zip(jobs, errors):
        if err:
            ledger.append(dict(base, file=f["name"], reason=err))
            continue
        blob = open(local, "rb").read()
        any_tree = False
        for member, data in files_in(f["name"], blob):
            m = measure_tree(data)
            if not m:
                continue
            any_tree = True
            measured.append((base, member, m))
        if not any_tree:
            ledger.append(dict(base, file=f["name"], reason="no tree read"))
    # classify every new genus once, in parallel, before screening
    need = {}
    for _, _, m in measured:
        if m["species"] and m["keyed_tips"] >= MIN_SPECIES_SHARE * m["tips"]:
            for k in m["species"]:
                need.setdefault(k.split("_")[0], []).append(k)
    todo = [(g, ex) for g, ex in need.items() if g not in clf.cache]
    log("classifying %d genera" % len(todo))
    with ThreadPoolExecutor(8) as pool:
        list(pool.map(lambda ge: clf.get(ge[0], ge[1]), todo))
    clf.save()
    for base, member, m in measured:
            s_ = screen(m, groups, basis, clf)
            if s_.get("verdict") and re.search(r"matri", member.split("!")[-1], re.I):
                # a Nexus data matrix can carry an analysis input tree (MrBayes
                # startvals): reported, not screened
                s_.update(verdict="", reason="tree inside a data matrix: an "
                          "analysis input, not a result")
            ledger.append(dict(base, file=member, tips=m["tips"],
                               keyed_tips=m["keyed_tips"],
                               species=len(m["species"]),
                               n_trees=m["n_trees"],
                               cv="" if m["cv"] == "" else round(m["cv"], 4),
                               ultrametric="" if m["cv"] == "" else
                               ("yes" if m["cv"] <= CHRONOGRAM_TOLERANCE
                                else "no"), intervals=m.get("intervals", ""), **s_))
    clf.save()
    ledger = compact(ledger)
    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(os.path.join(out_dir, "ledger.csv.gz"), "wt", newline="") as fh:
        w = csv.DictWriter(fh, LEDGER_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(ledger)
    with gzip.open(os.path.join(out_dir, "records.csv.gz"), "wt", newline="") as fh:
        cols = list(rec_rows[0].keys()) if rec_rows else ["source"]
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        w.writerows(rec_rows)
    write_report(since, until, records, rec_rows, ledger, out_dir)
    return ledger


def write_report(since, until, records, rec_rows, ledger, out_dir):
    by_src = Counter(r["source"] for r in records)
    with_trees = Counter(r["source"] for r in rec_rows if r["tree_files"])
    measured = [r for r in ledger if r.get("tips")]
    n_gene = sum(int(r.get("n_trees") or 0) for r in ledger
                 if str(r.get("reason", "")).endswith("trees not species-level"))
    species_level = [r for r in measured if r.get("group") or
                     r.get("reason") not in ("not species-level",)]
    cands = [r for r in measured if r.get("verdict")]
    blocked = [r for r in ledger if str(r.get("reason", "")).startswith("not downloadable")]
    title = ("# Living atlas search of the whole archive, %s" % until
             if since == "search" else "# Living atlas watch, %s to %s" % (since, until))
    L = [title, "",
         "Run %s. Rules and thresholds: `living/watch.py` docstring." %
         dt.datetime.now().strftime("%Y-%m-%d %H:%M"), ""]
    if HARVEST_ERRORS:
        L += ["**INCOMPLETE: the harvest missed records.** Rerun this window "
              "(`living/watch.py run --since %s --until %s`)." % (since, until), ""]
        L += ["- " + x for x in HARVEST_ERRORS] + [""]
    L += ["| Source | Records | With tree files |", "|---|---:|---:|"]
    for s in ("zenodo", "figshare", "dryad"):
        L.append("| %s | %d | %d |" % (s, by_src[s], with_trees[s]))
    best = OrderedDict()
    for r in sorted(cands, key=lambda r: -int(r["group_species"])):
        kinds = ["any"] + (["dated"] if "DATED" in r["verdict"] else [])
        for kind in kinds:
            best.setdefault((r["doi"], r["group"], kind), [r, 0])[1] += 1
    L += ["", "Trees measured: %d. Species-level: %d. Assigned to one group: %d. "
          "Candidate trees: %d, in %d deposits. Dryad tree files not "
          "downloadable: %d." % (
              len(measured) + n_gene,
              sum(1 for r in measured if r.get("reason") != "not species-level"),
              sum(1 for r in measured if r.get("group") and
                  not str(r.get("reason", "")).startswith("spans")),
              len(cands), len(best), len(blocked)), ""]
    L += ["## Candidates", ""]
    if cands:
        L += ["| Verdict | Group | Species in group | Not in reference | "
              "Secondary tree (accepted names) | Primary, dated (accepted names) | "
              "Tips | Trees | Ultrametric | Record | File |",
              "|---|---|---:|---:|---:|---:|---:|---:|---|---|---|"]
        shown = set()
        for (doi, grp, kind), (r, n) in sorted(best.items(),
                                               key=lambda kv: kv[0][1:] + kv[0][:1]):
            if id(r) in shown:
                continue
            shown.add(id(r))
            more = " (+%d more tree files)" % (n - 1) if n > 1 else ""
            L.append("| %s | %s | %s | %s | %s (%s) | %s (%s) | %s | %s | %s | [%s](https://doi.org/%s) %s | `%s` |" % (
                r["verdict"], r["group"], r["group_species"],
                r.get("not_in_reference", ""),
                r.get("reference") or "none", r.get("reference_species", 0),
                r.get("primary_reference") or "none", r.get("primary_species") or 0,
                r["tips"], r["n_trees"], r.get("ultrametric", "") +
                (" (intervals)" if r.get("intervals") and r.get("ultrametric") == "yes"
                 else ""),
                r["doi"], r["doi"], re.sub(r"<[^>]+>", "", r["title"])[:70]
                .replace("|", "/"), r["file"][-60:] + "`" + more + "`"))
    else:
        L.append("None in this window.")
    if blocked:
        L += ["", "## Dryad records with tree files (not measured)", "",
              "Dryad requires an API account for downloads. These records match the "
              "query and carry tree files:", ""]
        seen = OrderedDict()
        for r in blocked:
            seen.setdefault(r["doi"], (r["date"], r["title"], []))[2].append(r["file"])
        for doi, (date, title, files) in seen.items():
            L.append("- %s [%s](https://doi.org/%s) %s — %s" % (
                date, doi, doi, title[:90], ", ".join(files[:4])))
    L += ["", "## Next step", "",
          "Each candidate is curated before anything in the atlas changes: accepted "
          "names by `scripts/match_col.py`, off-group tips (D7), molecular tips, "
          "dating and archived date uncertainty from the source study, each with "
          "its evidence. A change is proposed on a branch and published only "
          "after approval."]
    with open(os.path.join(out_dir, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def _read_csv_any(base):
    for p, op in ((base + ".gz", gzip.open), (base, open)):
        if os.path.exists(p):
            with op(p, "rt", newline="") as fh:
                return list(csv.DictReader(fh)), p
    return [], None


def rebuild_report(run_dir):
    """Compact the ledger, rewrite it gzip-compressed, rebuild report.md."""
    since, until = os.path.basename(run_dir.rstrip("/")).split("_")
    ledger, lp = _read_csv_any(os.path.join(run_dir, "ledger.csv"))
    rec_rows, rp = _read_csv_any(os.path.join(run_dir, "records.csv"))
    ledger = compact(ledger)
    with gzip.open(os.path.join(run_dir, "ledger.csv.gz"), "wt", newline="") as fh:
        w = csv.DictWriter(fh, LEDGER_COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(ledger)
    if rp and not rp.endswith(".gz"):
        with gzip.open(rp + ".gz", "wt", newline="") as fh:
            w = csv.DictWriter(fh, list(rec_rows[0].keys()))
            w.writeheader()
            w.writerows(rec_rows)
        os.remove(rp)
    if lp and not lp.endswith(".gz"):
        os.remove(lp)
    records = [{"source": r["source"]} for r in rec_rows]
    write_report(since, until, records, rec_rows, ledger, run_dir)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--since")
    r.add_argument("--until", default=dt.date.today().isoformat())
    r.add_argument("--last-month", action="store_true")
    rp = sub.add_parser("report", help="rebuild a run's report from its files")
    rp.add_argument("run_dir")
    sp = sub.add_parser("search", help="search the whole archive, group by group, "
                                       "for trees not seen before")
    sp.add_argument("--groups", nargs="*", default=[],
                    help="group names (default: every eukaryotic group)")
    a = ap.parse_args(argv)
    if a.cmd == "report":
        rebuild_report(a.run_dir)
        return
    if a.cmd == "search":
        global SEARCH_MODE
        SEARCH_MODE = True
        today = dt.date.today().isoformat()
        out_dir = os.path.join(RUNS, "search_%s" % today)
        os.makedirs(out_dir, exist_ok=True)
        logf = open(os.path.join(out_dir, "run.log"), "a")

        def slog(msg):
            line = "%s %s" % (dt.datetime.now().strftime("%H:%M:%S"), msg)
            print(line, flush=True)
            logf.write(line + "\n")
            logf.flush()
        slog("search of the whole archive, groups: %s" % (" ".join(a.groups) or "all"))
        ledger = run("search", today, out_dir, slog,
                     harvest=search_harvest(set(a.groups)))
        slog("done: %d ledger rows, %d candidates; %s" % (
            len(ledger), sum(1 for x in ledger if x.get("verdict")),
            os.path.join(out_dir, "report.md")))
        if HARVEST_ERRORS:
            slog("INCOMPLETE search: %d problem(s); see report.md" % len(HARVEST_ERRORS))
            sys.exit(2)
        return
    if a.last_month:
        first = dt.date.today().replace(day=1)
        until = first - dt.timedelta(days=1)
        since, until = until.replace(day=1).isoformat(), until.isoformat()
    else:
        if not a.since:
            ap.error("--since or --last-month")
        since, until = a.since, a.until
    out_dir = os.path.join(RUNS, "%s_%s" % (since, until))
    os.makedirs(out_dir, exist_ok=True)
    logf = open(os.path.join(out_dir, "run.log"), "a")

    def log(msg):
        line = "%s %s" % (dt.datetime.now().strftime("%H:%M:%S"), msg)
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()
    log("window %s to %s" % (since, until))
    ledger = run(since, until, out_dir, log)
    log("done: %d ledger rows, %d candidates; %s" % (
        len(ledger), sum(1 for x in ledger if x.get("verdict")),
        os.path.join(out_dir, "report.md")))
    if HARVEST_ERRORS:
        log("INCOMPLETE harvest: %d problem(s); see report.md" % len(HARVEST_ERRORS))
        sys.exit(2)


if __name__ == "__main__":
    main()
