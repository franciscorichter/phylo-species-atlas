# Registers

The tables from which the paper's coverage levels, crown ages and sensitivity analyses are
computed. Each is written by a script in `code/rebuild/` or curated from the source publications
and cited there; `scripts/reproduce_paper.py` reads them, recomputes every reported number and
compares it with `paper/`.

## The per-group table

| File | What it holds |
|---|---|
| `coverage_basis.csv` | One row per group of described species (62: 60 eukaryotic, 2 prokaryotic). Described species and their authority, the reference tree, tips, accepted names (`c0_units`), molecular accepted names (`c1_units`) and how they were listed (`c1_basis`), the molecular fraction (`f_mol`), whether the tree is a chronogram in absolute time, and whether its source archive carries a counted crown interval. Groups without a reference tree are rows with zero tips. |
| `atlas_coverage.py`, `atlas_coverage.R` | A reader that forms any coverage figure as a numerator and denominator over `coverage_basis.csv`. |

## Reference trees

| File | What it holds |
|---|---|
| `canonical_sources.csv` | For each reference tree, the unpruned source it is built from (`standardized/canonical_sources/`) and where that source came from. |
| `canonical_source_prep.csv`, `canonical_source_crosswalks/` | How the sources taken from a deposit other than an earlier redistributed file were relabelled: deposited file, tree index, branch-length scale, label crosswalk. |
| `pruned_tips.csv` | Tips outside their group (outgroups, misplaced taxa), removed before counting, with the reason for each. |
| `offtarget_checks.csv` | For each reference tree, how off-group tips were searched for and what was found. |
| `tree_provenance.csv` | Per source tree: citation, DOI, data type, inference, dating method, Open Tree pipeline involvement, imputation status. |
| `dataset_table_registry.csv` | The rows of the paper's dataset table: citation key, year, journal and flags for every redistributed tree. |

## Molecular tips

| File | What it holds |
|---|---|
| `molecular_tips.csv` | Tip ids of each reference tree that are in the molecular analysis behind it, for the 23 groups (21 eukaryotic, Archaea and Bacteria) whose molecular tips are listed. In the other three groups every tip is in the molecular analysis, as their sources state. |
| `molecular_sources/<group>.taxa.txt` | The taxa with at least one site of sequence data in a source study's alignment (13 groups). The header names the alignment, its DOI or URL and its checksum. Written by `code/rebuild/alignment_taxa.py`. |
| `molecular_sources/source_trees/<group>.taxa.txt` | The tip labels of the tree of sequenced species that a source study distributes (7 groups). The trees themselves are not redistributed; the header names the file, the study and the checksum. Written by `code/rebuild/export_source_tree_taxa.py`. |
| `molecular_tip_evidence.csv` | Tip-by-tip evidence of sequence data for the three groups whose sources give no machine-readable list (Cnidaria, Crustaceans, Echinoderms). |
| `molecular_only_trees.csv` | For each group whose reference tree holds taxonomically placed tips, the tree of molecular species its source distributes. |

## Dating and intervals

| File | What it holds |
|---|---|
| `time_scale_checks.csv` | For each ultrametric reference tree, whether its branch lengths are in absolute or relative time, and the evidence. |
| `dating_deposit_checks.csv` | What each source study dated and what its archive holds (`study_dating`). |
| `crown_clades.csv` | The focal clade whose crown age is reported for each chronogram. |
| `crown_interval_sources.csv` | The archived tree set or annotated tree from which each crown interval is read, its kind and burn-in. |
| `crown_intervals.csv` | The crown interval of each source archive, its kind, and the tests that decide whether it counts (degenerate, at a calibration bound, reproducing the calibration prior). Written by `code/rebuild/crown_intervals.py`. |
| `calibration_register.csv` | The calibration each source placed on the focal node, with its prior range and the quotation it is taken from. |

## Denominators and name authorities

| File | What it holds |
|---|---|
| `col_denominators.csv` | Catalogue of Life (COL26.9 base release, and the extended release XR) species counts per group, the single-authority denominator of the sensitivity analysis. Written by `code/rebuild/col_denominators.py`. |
| `external/gbif_unmatched.csv` | GBIF backbone matches of the tip names the Catalogue of Life leaves unmatched (`code/rebuild/gbif_unmatched.py`), with the query date. |
| `external/opentree16.1_counts.csv`, `opentree16.1_phylo_species.tsv.gz`, `opentree16.1_phylo_species_col.csv` | Open Tree of Life synthesis 16.1: tips and phylogeny-placed species per group, the species list, and its Catalogue of Life resolution (`code/rebuild/opentree_extract.py`, `match_col.py --opentree`). |
| `opentree_clades.csv` | The Open Tree taxonomy ids of each group's clade. |

## Inputs of the rebuild

`_rebuild/truth_*.csv` are the per-file and per-group fields that `code/rebuild/rebuild_atlas.py`
reads: group of every redistributed file, measured tip counts and root-to-tip variation, described
species and estimate ranges, and the provenance fields of the reference trees.
