# Standardized trees and names

The redistributed trees, with numeric tip ids, and the name layer that resolves the ids.

| File | What it holds |
|---|---|
| `trees/` | One Newick file per redistributed tree; tip labels are numeric ids. The 26 reference trees (one per group with a tree) are built by rule from `canonical_sources/` (below); every other file is a copy of a source deposit with its tips relabelled. |
| `metadata.csv` | Per tree: file name, group, study, number of tips (counted in the file), and the `dated` flag of the earlier releases. Whether a tree is a chronogram is measured from the file instead: `registers/_rebuild/truth_dated.csv` and `scripts/reproduce_paper.py`. |
| `dictionary.csv` | Tip id -> standardized name. It holds exactly the labels of the files in `trees/`, `canonical_sources/` and `superseded_canonicals/`. Ids are stable across releases: an id that occurs in this and an earlier release names the same label. |
| `full_mapping.csv` | Tip id, label as written in the source deposit, standardized name, group. The bridge from a source's labels to ids. |
| `col_resolution.csv.gz` | For every tip of every eukaryotic reference tree: the query formed from its name (the species key or a documented override), and its match in the Catalogue of Life COL26.9 base release (ChecklistBank 316321, doi:10.48580/dgz5n): match type, accepted-name id and name, kingdom, family. Named coverage counts the distinct accepted names whose match type resolves (`accepted`, `provisional`, `synonym` and their orthographic variants). |
| `col_resolution_noncanonical.csv` | The same resolution for the tips of every other redistributed eukaryotic tree (the union sensitivity). |
| `col269_family_species_counts.csv` | Accepted species per family in COL26.9 (the family sampling fractions of the worked example). |
| `col269_suprageneric_genera.txt` | Genera of COL26.9 whose names end in a family-level suffix (`Xerophyta`, `Elapsoidea`, ...), read by the species-key rule so that these are not taken for family names. |
| `name_overrides.csv` | Tips whose label carries no binomial (lineage codes) or a documented misspelling, with the query used and its evidence. |
| `canonical_sources/` | The unpruned, id-labelled source of each reference tree. |
| `superseded_canonicals/` | The former reference trees of seven groups (Birds, Bryophytes, Diatoms, Fish, Gastropods, Sharks, Turtles), replaced by a tree the same source deposited: a chronogram where the study archived one, the tree with usable branch lengths (Diatoms), the full all-taxon tree (Fish), the maximum-clade-credibility tree of the posterior set (Sharks). Kept so that every label of the dictionary occurs in a distributed file. |
| `tree_references.csv` | For every redistributed tree and canonical source, the file its branch lengths were checked against and the rule of the check. |
| `withdrawn.csv` | Files not redistributed (the TimeTree of Life, at its project's request), with checksums of every known copy. |

## How a reference tree is built

From its source in `canonical_sources/` (`registers/canonical_sources.csv`):

1. the tips listed in `registers/pruned_tips.csv` (outside the group) are removed;
2. walking the tips left to right, a tip whose species key was already seen is removed, so one
   tip per species remains; the species key is the first `Genus_epithet` pair of the
   standardized name, and a name with cf., aff., sp. or a hybrid marker has no key and is kept
   as an unresolved label;
3. nodes left with one child are collapsed, their branch lengths summed, so the root is the most
   recent common ancestor of the retained tips.

`scripts/reproduce_paper.py` checks, for every reference tree, that its tips are exactly those
the rule keeps and that every root-to-tip depth equals the source's less one constant.
`code/rebuild/canonical_trees.py` is the script that writes the files.

## Loading a tree

```r
library(ape)
tree <- read.tree("standardized/trees/mammals.nwk")
dict <- read.csv("standardized/dictionary.csv")
tree$tip.label <- dict$standardized_name[match(as.integer(tree$tip.label), dict$id)]
```
