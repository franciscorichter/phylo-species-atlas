#!/usr/bin/env python3
"""alignment_taxa.py -- list the taxa with sequence data in a deposited alignment.

Writes the files of registers/molecular_sources/<group>.taxa.txt: the taxa of a
source study's alignment that carry at least one site that is neither a gap
nor missing data.  Standard library only.

    python3 alignment_taxa.py FORMAT INPUT [options] > group.taxa.txt

FORMAT       fasta | phylip | nexus | clustal | beast-xml | iqtree-log | csv
INPUT        the alignment as deposited (.gz is read directly); for a file
             inside a zip archive give ARCHIVE.zip --member NAME
--alphabet   dna (gap and missing: - ? . N n) or protein (- ? . X x *);
             default dna
--expand-prefixes TREE
             CLUSTAL truncates names; a name that is not a tip label of TREE
             is replaced by the one tip label of TREE it is a prefix of
--column C --value V --label-column L
             csv: list column L of the rows where column C equals V
--header TEXT
             a comment line written first (repeatable): source, DOI, checksum

Formats as the sources deposit them:

    fasta       Bivalves (Pfeiffer et al. 2019), Diatoms (Nakov et al. 2018),
                Ferns (FTOL v1.8.0), Insects (insectphylo Insecta_v2),
                Archaea and Bacteria (GTDB R10-RS226, --alphabet protein)
    phylip      Bryophytes (Bechteler et al. 2023), Gastropods (Zapata et al.
                2014, --alphabet protein)
    nexus       Bryozoa (Orr et al. 2022, TRIMMED.nex in the Dryad zip)
    beast-xml   Cephalopods (Basava 2024): a taxon is listed when any locus
                has a site with data
    clustal     Sponges (Lavrov et al. 2023, --expand-prefixes)
    iqtree-log  Fungi (Li et al. 2021): the per-sequence table of the IQ-TREE
                log; a taxon is listed when its gap/ambiguity share is below
                100%
    csv         Mammals (Upham et al. 2019): taxonomy_mamPhy_5911species_
                toPublish.csv, --column samp --value sampled
                --label-column tiplabel
"""

import argparse
import csv
import gzip
import io
import re
import sys
import zipfile

MISSING = {"dna": set("-?.Nn"), "protein": set("-?.Xx*")}


def open_text(path, member=None):
    if member:
        with zipfile.ZipFile(path) as z:
            return io.TextIOWrapper(io.BytesIO(z.read(member)), encoding="utf-8",
                                    errors="replace")
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def has_data(seq, missing):
    return any(ch not in missing and not ch.isspace() for ch in seq)


def fasta(fh, missing):
    name, data, out = None, False, []
    for ln in fh:
        ln = ln.rstrip("\n")
        if ln.startswith(">"):
            if name is not None and data:
                out.append(name)
            name, data = ln[1:].strip().split()[0] if ln[1:].strip() else "", False
        elif name is not None and not data and has_data(ln, missing):
            data = True
    if name is not None and data:
        out.append(name)
    return out


def phylip(fh, missing):
    lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    ntax = int(lines[0].split()[0])
    body = lines[1:]
    names, seqs = [], []
    for ln in body[:ntax]:
        parts = ln.split(None, 1)
        names.append(parts[0])
        seqs.append(parts[1] if len(parts) > 1 else "")
    for i, ln in enumerate(body[ntax:]):           # interleaved blocks
        seqs[i % ntax] += ln
    return [n for n, s in zip(names, seqs) if has_data(s, missing)]


def nexus(fh, missing):
    seqs, order, inmat = {}, [], False
    for ln in fh:
        s = re.sub(r"\[[^\]]*\]", "", ln).strip()
        if not inmat:
            if s.lower().startswith("matrix"):
                inmat = True
            continue
        if s.startswith(";"):
            break
        if not s:
            continue
        end = s.endswith(";")
        s = s.rstrip(";").strip()
        parts = s.split(None, 1)
        if parts:
            nm = parts[0].strip("'")
            if nm not in seqs:
                seqs[nm] = ""
                order.append(nm)
            seqs[nm] += parts[1] if len(parts) > 1 else ""
        if end:
            break
    return [n for n in order if has_data(seqs[n], missing)]


def clustal(fh, missing):
    seqs, order = {}, []
    for ln in fh:
        if ln.startswith("CLUSTAL") or not ln.strip() or ln[0].isspace():
            continue
        parts = ln.split()
        if len(parts) < 2:
            continue
        if parts[0] not in seqs:
            seqs[parts[0]] = ""
            order.append(parts[0])
        seqs[parts[0]] += parts[1]
    return [n for n in order if has_data(seqs[n], missing)]


def beast_xml(fh, missing):
    found, order = {}, []
    rx = re.compile(r'<sequence\b[^>]*\btaxon="([^"]+)"[^>]*\bvalue="([^"]*)"')
    for m in rx.finditer(fh.read()):
        t, v = m.group(1), m.group(2)
        if t not in found:
            found[t] = False
            order.append(t)
        found[t] = found[t] or has_data(v, missing)
    return [t for t in order if found[t]]


def iqtree_log(fh, missing):
    out, inside = [], False
    rx = re.compile(r"^\s*\d+\s+(\S+)\s+([\d.]+)%")
    for ln in fh:
        if "Gap/Ambiguity" in ln:
            if out:
                break                       # the first table only
            inside = True
            continue
        if inside:
            m = rx.match(ln)
            if m:
                if float(m.group(2)) < 100.0:
                    out.append(m.group(1))
            elif out:
                break
    return out


def newick_labels(path):
    with open_text(path) as fh:
        text = fh.read()
    text = re.sub(r"\[[^\]]*\]", "", text)
    return set(t.strip().strip("'") for t in
               re.findall(r"[(,]\s*([^(),:;\s]+)", text))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("format", choices=["fasta", "phylip", "nexus", "clustal",
                                       "beast-xml", "iqtree-log", "csv"])
    ap.add_argument("input")
    ap.add_argument("--member")
    ap.add_argument("--alphabet", choices=["dna", "protein"], default="dna")
    ap.add_argument("--expand-prefixes")
    ap.add_argument("--column")
    ap.add_argument("--value")
    ap.add_argument("--label-column")
    ap.add_argument("--header", action="append", default=[])
    a = ap.parse_args(argv)
    missing = MISSING[a.alphabet]
    with open_text(a.input, a.member) as fh:
        if a.format == "csv":
            taxa = [r[a.label_column] for r in csv.DictReader(fh)
                    if r[a.column] == a.value]
        else:
            taxa = {"fasta": fasta, "phylip": phylip, "nexus": nexus,
                    "clustal": clustal, "beast-xml": beast_xml,
                    "iqtree-log": iqtree_log}[a.format](fh, missing)
    if a.expand_prefixes:
        labels = newick_labels(a.expand_prefixes)
        out = []
        for t in taxa:
            if t in labels:
                out.append(t)
                continue
            hits = sorted(l for l in labels if l.startswith(t))
            if len(hits) != 1:
                sys.exit("%s: %d tree labels start with it; cannot expand"
                         % (t, len(hits)))
            out.append(hits[0])
        taxa = out
    w = sys.stdout
    for h in a.header:
        w.write("# %s\n" % h)
    for t in sorted(set(taxa)):
        w.write(t + "\n")


if __name__ == "__main__":
    main()
