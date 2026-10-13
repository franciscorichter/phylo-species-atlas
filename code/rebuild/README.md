# The pipeline that produced the release

These scripts wrote the files of `standardized/`, `registers/` and `paper/`. They are copied
unchanged from the paper's working repository, where the files of this repository sit under
`Data/` and `Paper/inputs/`, and the source deposits the studies archived sit under `Data/raw/`.
The source deposits are not redistributed: they remain in their original archives, cited in
`registers/tree_provenance.csv` and in the headers of `registers/molecular_sources/`.

`release/path_map.csv` of the working repository lists every file of this repository that the
pipeline writes, with its path there; `scripts/build_release.py` copies them here.

`scripts/reproduce_paper.py` is independent of these scripts: it recomputes the paper's numbers
from this repository by the stated rules and does not import them.

## Path map

| In the scripts | In this repository |
|---|---|
| `Data/standardized/trees/`, `dictionary.csv`, `full_mapping.csv`, `metadata.csv`, ... | `standardized/` (same names) |
| `Data/standardized/col_resolution.csv` | `standardized/col_resolution.csv.gz` (gzip) |
| `Data/standardized/names/<file>.json` | not stored; `site/build.py` writes the same id -> name maps from `dictionary.csv` |
| `Data/_archive/2026-10-revision/canonical_sources/` | `standardized/canonical_sources/` |
| `Data/_archive/2026-10-revision/superseded_canonicals/trees/` | `standardized/superseded_canonicals/` |
| `Data/<register>.csv`, `Data/molecular_sources/`, `Data/canonical_source_crosswalks/`, `Data/external/`, `Data/_rebuild/` | `registers/` (same names) |
| `Paper/inputs/` | `paper/` |
| `Paper/audit-2026-08-20/code/common.py` | `code/rebuild/audit_common.py` |
| `Paper/code/analysis/worked_example_bd.R` | `code/rebuild/worked_example_bd.R` |
| `Data/_archive/2026-10-revision/corrupted_branch_lengths/RESTORED.csv` | `registers/restored_branch_lengths.csv` |
| `Data/raw/` | not redistributed |

## Scripts

| Script | Writes | Needs |
|---|---|---|
| `canonical_trees.py` | the reference trees in `standardized/trees/` from `standardized/canonical_sources/` by the rule in `standardized/README.md` | this repository |
| `prepare_canonical_sources.py` | the canonical sources taken from a deposit (`registers/canonical_source_prep.csv`) | the deposits |
| `restore_branch_lengths.py` | the trees with restored branch lengths (`ERRATUM.md`; first released in v1.1.1) | the pre-relabelling copies |
| `alignment_taxa.py` | `registers/molecular_sources/<group>.taxa.txt` from the cited alignments; its docstring gives the call for each group | the alignments (GTDB, Dryad, Zenodo, Figshare, GitHub, as cited in each header) |
| `export_source_tree_taxa.py` | `registers/molecular_sources/source_trees/` from the trees of sequenced species | the source trees |
| `match_col.py` | `standardized/col_resolution*.csv`, `col269_*` | the COL26.9 release (ChecklistBank 316321) |
| `col_denominators.py` | `registers/col_denominators.csv` | the ChecklistBank API |
| `gbif_unmatched.py` | `registers/external/gbif_unmatched.csv` | the GBIF API |
| `opentree_extract.py` | `registers/external/opentree16.1_*` | Open Tree synthesis 16.1 and OTT 3.7.3 |
| `crown_intervals.py` | `registers/crown_intervals.csv` | the archived tree sets and annotated trees (`registers/crown_interval_sources.csv`) |
| `summary_tree_node_ages.py` | `registers/summary_tree_node_ages.csv` | Python 3; the AvesData 1.3 summary tree, `all_node_ages.json` and `phylo_only_clements_labels.tre`, not redistributed |
| `taxon_addition_dates.R` | `registers/taxon_addition_set_check.csv` | R with ape; the fish tree set and dated tree of FishTreeOfLife, not redistributed |
| `worked_example_bd.R` | `paper/worked_example_bd*.csv` except `worked_example_bd_sets.csv` | R with ape (phytools for the cross-check); this repository. `run_worked_example.py` runs it from this repository and compares the result |
| `worked_example_sets.R` | `paper/worked_example_bd_sets.csv` | R with ape; the archived tree sets of `registers/worked_example_tree_sets.csv` (VertLife, AvesData), not redistributed; see `README.md` |
| `rebuild_atlas.py` | `registers/coverage_basis.csv`, `molecular_tips.csv` and every file of `paper/` | the paper repository layout, PyYAML |
| `verify_atlas.py` | nothing; the acceptance test of the paper repository (re-derives every generated file) | the paper repository layout |

Standard library Python unless stated; `rebuild_atlas.py` needs PyYAML.
