# Erratum: branch lengths in releases v1.0.3 to v1.0.8

**Affected releases:** v1.0.3 (2026-05-12) through v1.0.8 (2026-07-05), on GitHub and on Zenodo
(concept DOI 10.5281/zenodo.20127157). **Corrected in:** v1.1.1 (10.5281/zenodo.23455934) and every later release.
Releases v1.0.0 to v1.0.2 are not affected.

The Zenodo record 10.5281/zenodo.23348695, labelled v1.1.0 with publication date 2026-10-13,
holds the v1.0.8 files, that is the corrupted branch lengths; it is not a release of this
repository's v1.1.x line, and its metadata says so. There is no release v1.1.0: the first
release with corrected branch lengths is v1.1.1.

## What was wrong

The step that relabelled tips with numeric dictionary ids, introduced in v1.0.3, rewrote digit
runs elsewhere in the Newick strings. A digit run immediately followed by `,`, `)` or `:` was
treated as a tip id: it was replaced through the old-to-new id map, or it lost its leading zeros.
This changed:

- **branch lengths**, e.g. `0.007143` became `0.7146` and `0.015603` became `0.15695`;
- **numeric internal-node labels** followed by `:`, including support values and stored age
  intervals, e.g. `11.1189_2.0251` became `11.1189_2.251`.

Tip sets and topologies are unchanged. Root-edge lengths (followed by `;`) are unchanged.

Every redistributed tree that carried the defect is corrected in v1.1.1. The register
`registers/restored_branch_lengths.csv` lists every file the restoration checked, one row per file,
with `changed` = `yes` where the file differed from its pre-relabelling copy. Its `file` column
gives the path in the repository that builds this release: `Data/standardized/trees/<name>` is
`standardized/trees/<name>` here, and `Data/site/...` is `site/...`. Of its 306 rows, 299 have
`changed` = `yes`: 256 are trees redistributed in `standardized/trees/`, 1 is `timetree.nwk`,
which was restored and is no longer redistributed (`standardized/withdrawn.csv`), and 42 are the
copies of trees under `site/data/partitions/`, of which the browser site holds 18. The site holds
no copy of a reference tree: it reads every redistributed tree from `standardized/trees/`, and its
Primates tree with one tip per species is derived from the restored `standardized/trees/primates.nwk`.
The other 7 rows are files of `standardized/trees/`
that were checked and found unchanged. Any analysis that used branch lengths, node ages or support values from these files in
an affected release should be repeated with v1.1.1 or later. Analyses that used only tip sets or
topologies are not affected.

## How it was corrected

Each affected file was rebuilt from its pre-relabelling copy, replacing tip ids only and copying
every other token verbatim. Every corrected file was checked to have the same tip multiset and the
same topology as the file it replaces. The corrected branch lengths were checked against the
original source files deposited by the studies that published the trees.

## Effect on reported results

The first version of the preprint (bioRxiv doi:10.64898/2026.05.31.729127, posted 31 May 2026)
was written while releases v1.0.3 onward were current; read its branch-length results, and any
files attached to it, against v1.1.1.

Measures computed from branch lengths change: root-to-tip variation, the classification of trees
as ultrametric, and crown ages. Counts of species and tips do not. Revised figures are reported in
the next version of the accompanying manuscript.
