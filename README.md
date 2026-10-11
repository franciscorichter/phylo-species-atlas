# Phylogenetic Species Atlas

Published species-level phylogenies of eukaryotes, redistributed with numeric tip ids, one
name layer resolved to the Catalogue of Life, and the registers from which the coverage of
described species by these trees is counted. It is the data and code of the paper

> Richter F. *Described eukaryotic species in published species-level phylogenies: placement by molecular data, dates and date uncertainty.* bioRxiv,
> doi:[10.64898/2026.05.31.729127](https://doi.org/10.64898/2026.05.31.729127).

Every number in the text, tables and figures of the paper and of its Supplement is recomputed
from this repository by

```sh
python3 scripts/reproduce_paper.py
```

> **Branch lengths in releases v1.0.3 to v1.0.8 are corrupted.** Use v1.1.0 or later for any
> analysis of branch lengths, node ages or support values. Details: [`ERRATUM.md`](ERRATUM.md).
>
> Two files carry the literal value `NaN` as a branch length: `standardized/trees/ants.nwk`
> (83 branches) and `standardized/superseded_canonicals/diatoms.nwk` (95 branches). Some Newick
> readers reject it; neither file is a reference tree.

Zenodo concept DOI (all versions): [10.5281/zenodo.20127157](https://doi.org/10.5281/zenodo.20127157).
Licence: CC BY 4.0 (`LICENSE`); the source trees keep the terms of their original archives.

## Contents

| Path | What it holds |
|---|---|
| `standardized/trees/` | 263 Newick files with numeric tip ids: the 26 reference trees (one per group that has a species-level tree, 24 eukaryotic and 2 prokaryotic) and 237 other trees, among them 218 family-level trees of tetrapods (Condamine et al. 2019), one of which is the reference tree of Crocodilians. |
| `standardized/` | The name layer: `dictionary.csv` (id -> standardized name), `full_mapping.csv` (source label -> name), the Catalogue of Life resolution of every reference-tree tip (`col_resolution.csv.gz`), per-tree `metadata.csv`, `canonical_sources/` (the unpruned source of each reference tree) and `superseded_canonicals/` (the seven reference trees this release replaces). See [`standardized/README.md`](standardized/README.md). |
| `registers/` | The per-group table (`coverage_basis.csv`), the molecular tips of each reference tree and the lists they are taken from, the time-scale, crown-interval and calibration registers, the pruned off-group tips, and the Catalogue of Life, GBIF and Open Tree extracts. See [`registers/README.md`](registers/README.md). |
| `paper/` | The paper's generated numbers, tables and figure data, as the manuscript prints them. `paper/revision/` holds the review register of the revision (`FINDINGS.tsv`) and its recorded decisions (`DECISIONS.md`), described in the paper's Use of AI section. |
| `scripts/reproduce_paper.py` | Recomputes the paper's numbers from `standardized/` and `registers/` and compares them with `paper/`. Python 3, standard library only. |
| `code/rebuild/` | The pipeline that wrote `standardized/`, `registers/` and `paper/`, with a map from its paths to this repository. See [`code/rebuild/README.md`](code/rebuild/README.md). |
| `code/figures/` | The scripts of the five main-text figures (matplotlib). |
| `phyloatlas/` | An R package that loads any tree by name. |
| `site/`, `audits/`, `data_provenance.csv`, `data_estimates.csv`, `data_sources.csv` | The browser site and its inputs (see *The browser site* below). |

## The coverage levels

The 2,085,869 described eukaryotic species are divided into 60 groups; 24 of them have a
species-level phylogeny, and for each the largest one found is the group's reference tree, except
Birds: AvesData 1.3 aligns one synthesis to three versions of the Clements taxonomy, and the
version used is the one with the most tips placed by an input phylogeny (Clements 2021, 9,239
such tips; the larger Clements 2023 version has 11,017 tips, 9,189 of them placed by an input
phylogeny, and is listed in `registers/larger_trees_not_counted.csv`). A
described species is counted at four nested levels:

| Level | A species counts when | Species | % of described |
|---|---|---:|---:|
| Named | its Catalogue of Life (COL26.9) accepted name is a tip of its group's reference tree | 395,517 | 19.0 |
| Molecular | it is named, and the tip is in the molecular analysis from which the tree was built | 157,201 | 7.5 |
| Dated | it is molecular, and the reference tree is a chronogram in absolute time | 113,180 | 5.4 |
| Uncertainty available | it is dated, and the source archive carries a node-age distribution or interval a user can propagate (a tree set whose node ages vary, or intervals on a summary tree), whatever its relation to the calibration prior | 23,942 | 1.1 |

A species named in two groups' trees counts once. Bacteria and Archaea are counted separately, in GTDB
species clusters, and enter no eukaryotic number. The values above are those of
`paper/headline_numbers.json`, which `scripts/reproduce_paper.py` recomputes.

How each level is decided, tip by tip:

- **Named**: `standardized/col_resolution.csv.gz` gives each tip's Catalogue of Life match; tips
  resolving to one accepted name count once; tips labelled cf., aff., sp. or as hybrids are
  unresolved labels.
- **Molecular**: `registers/molecular_tips.csv` lists the molecular tips of 21 eukaryotic
  reference trees, re-derivable from the alignment taxa lists and source-tree tip lists in
  `registers/molecular_sources/` and the per-tip evidence in `registers/molecular_tip_evidence.csv`;
  in the other three groups every tip is in the molecular analysis.
- **Dated**: the reference tree's root-to-tip coefficient of variation is at most 0.05, measured
  on the file, and `registers/time_scale_checks.csv` records its branch lengths as absolute time
  (the insect tree is ultrametric in relative time and does not count).
- **Uncertainty available**: `registers/crown_intervals.csv` (`uncertainty_available`) records
  whether the source archive carries a posterior, bootstrap, date-assignment or root-prior tree
  set, or node intervals on a summary tree. The fish taxon-addition set does not count: pruned to
  its sequenced tips, each of its trees has the dates of the source's one dated tree
  (`registers/taxon_addition_set_check.csv`). A stricter reading, printed beside it, counts only
  archives holding dated trees of a posterior or bootstrap dating of sequence data whose crown age
  is neither fixed nor on a bound (`kind`, `degenerate`, `at_bound`).

`registers/crown_intervals.csv` also compares each archived crown interval with the range its
calibrations allow (`not_excluded_by_rule`, `rule_exclusion`; `registers/calibration_register.csv`).
That comparison is a description, printed in Supplement S6, and no level rests on it: how far
the data moved an interval needs a run of the source's analysis that samples from the prior,
which was made for one group (Cephalopods).

## Reproducing the paper

```sh
python3 scripts/reproduce_paper.py            # under a minute; exits 0 when every number agrees
python3 scripts/reproduce_paper.py --quiet    # failures and summary only
```

The script re-implements the counting rules from their statements, without importing the
pipeline, and checks:

- every reference tree against its source: tips = source tips minus off-group tips minus repeated
  species, branch lengths preserved;
- the molecular tips of all 21 listed eukaryotic groups and both prokaryotic groups, re-derived
  from the distributed lists;
- the per-group table, the four levels (per group, per category, on the union of groups), the
  losses between levels, every sensitivity analysis the paper prints, the crown ages and the interval
  classification, the Open Tree comparison, and the inputs of the worked example;
- every number printed in the text (all macros of `paper/*_macros.tex`), the three main-text
  tables cell by cell, and the data of Figures 1 to 5;
- every table of Supplement S1 to S10 cell by cell (`paper/supp_*.csv`, with the counts its
  `.tex` twin prints) and every number the Supplement's text prints (`paper/supp_macros.tex`).

The birth-death fits of the worked example are made in R:

```sh
python3 code/rebuild/run_worked_example.py    # R with ape (and phytools); a few minutes
```

It runs `code/rebuild/worked_example_bd.R` on this repository's files and then checks every
printed worked-example number against the fresh fits.

The fits repeated over the dated tree sets the sources archive (`paper/worked_example_bd_sets.csv`,
Supplement S9 and Table S13) need those sets, which are not redistributed here, so
`run_worked_example.py` does not rerun them and checks that block against the committed copy.
To refit them, download the sets listed in `registers/worked_example_tree_sets.csv` into the paths
its `source` column gives under a `Data/` directory beside a copy of this repository's
`registers/` and `standardized/` files, and run `code/rebuild/worked_example_sets.R` (R with ape):

| Group | Set | Where it is archived |
|---|---|---|
| Amphibians | posterior, 10,000 trees | VertLife, `data.vertlife.org/amphibiantree` (Jetz & Pyron 2018) |
| Squamates | posterior, 10,000 trees | VertLife, `data.vertlife.org/squamatetree` (Tonini et al. 2016) |
| Mammals | NDexp completed and DNA-only sets, 10,000 trees each | VertLife, `data.vertlife.org/mammaltree` (Upham et al. 2019) |
| Birds | date-assignment sets, 100 trees each | AvesData release 1.3, `github.com/McTavishLab/AvesData`, `Tree_versions/Aves_1.3/Clements2021` (McTavish et al. 2025) |
| Sharks | 10-fossil set (10,000 trees) and the 500-tree set of sequenced species of the same scenario | VertLife, `data.vertlife.org/sharktree` (Stein et al. 2018) |

The figures: `python3 code/figures/fig_coverage_cascade.py` and the four other scripts in
`code/figures/` write Figures 1 to 5 to `paper/figures/`.

## What is not redistributed

- The source deposits (Dryad, Zenodo, Figshare, GitHub, VertLife, FishTreeOfLife, insectphylo,
  GTDB). Each is cited in `registers/tree_provenance.csv`, in the headers of
  `registers/molecular_sources/` and in `registers/crown_interval_sources.csv`. What the paper
  needs from them is redistributed as derived files: the relabelled trees, the taxa lists of the
  alignments, the tip lists of the trees of sequenced species, and the crown intervals read from
  the archived tree sets. `code/rebuild/` holds the scripts that derive them.
- The TimeTree of Life (Hedges et al. 2015), at the TimeTree project's request; it is available
  from timetree.org. `standardized/withdrawn.csv` lists the checksums of its copies.

## Loading a tree

With ape, from a clone of this release:

```r
library(ape)
tree <- read.tree("standardized/trees/mammals.nwk")
dict <- read.csv("standardized/dictionary.csv")
tree$tip.label <- dict$standardized_name[match(as.integer(tree$tip.label), dict$id)]
```

The `phyloatlas` R package loads a tree by name. It reads the files under the location named
by the option `phyloatlas.base_url`, which defaults to the `main` branch of this repository on
GitHub, not to this release. To read the files of this release, set the option to the root of a
clone of it (a directory path; the package reads `standardized/` and `data_provenance.csv`
under it) before loading:

```r
pak::pkg_install("franciscorichter/phylo-species-atlas/phyloatlas")
library(phyloatlas)
options(phyloatlas.base_url = "path/to/phylo-species-atlas")  # a clone of this release
list_trees()
tree <- load_atlas_tree("mammals")
```

The tag `v1.1.0` is created on GitHub when this release is published. From then on the option
can name the tag instead of a clone,
`"https://raw.githubusercontent.com/franciscorichter/phylo-species-atlas/v1.1.0"`; before
publication that address returns HTTP 404.

## The browser site

`site/` is the source of the page at
<https://franciscorichter.github.io/phylo-species-atlas/>, built by `site/build.py` from
`standardized/metadata.csv`, `data_provenance.csv`, `data_estimates.csv`, `data_sources.csv`,
`audits/` and the per-group notes in `site/data/partitions/`. Its content does not follow the
paper:

- its coverage percentages are tips over described species of each dataset, taken from
  `data_provenance.csv`, not the four levels above, and its described-species totals and crown
  ages are the values typed there;
- built from this release, it shows and offers for download the copies in
  `site/data/partitions/*/tree.nwk` instead of the files in `standardized/trees/`. Of the 24
  copies that share a name with a release file, all of them reference trees, 4 (Archaea,
  Bacteria, Bivalves, Cephalopods) are identical to it, 4 (Birds, Bryophytes, Sharks, Sponges)
  have its tips but differ in branch lengths or labels, and 16 have a different set of tips; the
  Cnidaria copy has no release file of that name (the Cnidaria reference tree is `corals.nwk`),
  and the Crocodilians reference tree has no copy;
- the committed `site/data.json` was built from release v1.0.8 and lists `timetree.nwk`, which is
  not redistributed. The deployment workflow (`.github/workflows/pages.yml`) rebuilds it from the
  files above on every push to `main`.

The coverage of described species by these trees is the one in `registers/coverage_basis.csv` and
`paper/`.
