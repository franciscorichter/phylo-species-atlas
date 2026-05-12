# Standardized archive

This directory contains the standardized tree archive described by the manuscript supplement. It is the publication-facing archive layer of the `phylo_trees` workspace.

---

## Intended contract

The standardized archive is supposed to provide a stable, compact representation of the representative phylogenetic datasets used by the project.

Core components:

- `trees/`
  - standardized Newick trees with numeric tip labels
- `dictionary.csv`
  - global mapping from numeric ID to standardized label
- `metadata.csv`
  - per-tree metadata
- `full_mapping.csv`
  - mapping from original labels to standardized labels and numeric IDs

According to the current manuscript supplement, the archive should contain:

- 264 standardized Newick trees in the current working archive
- 637,619 dictionary entries
- 702,787 full-mapping entries

These counts reflect the current reconstructed archive state. Older manuscript-era counts should be treated as historical audit targets rather than present-tense assumptions.

---

## Current state

The standardized archive has now been re-audited after reconstructing the previously missing representative groups `sharks`, `bees`, `salamanders`, `cetaceans`, and `bivalves`.

Current high-level status:

- all previously reconstructed representative gaps are filled, including `bryozoa` and `frogs` which were added in the latest audit cycle
- no remaining `standardized_missing` provenance groups
- file presence, metadata coverage, and mapping coverage are synchronized within the current working archive

Remaining work is mostly cleanup rather than gap-filling:

- move legacy figure byproducts out of this folder or classify them more explicitly as byproducts
- review the remaining alias and tip-count mismatch cases surfaced in the audit crosswalks
- keep manuscript-facing reported counts aligned with the reconstructed archive

---

## File roles

### `trees/`

Contains the standardized Newick trees that are intended to be shared publicly.

Properties:

- one file per representative dataset or Condamine family tree
- tip labels are numeric IDs, not species names
- the numeric IDs resolve through `dictionary.csv`

### `dictionary.csv`

Global two-column table:

- `id`
- `standardized_name`

This is the canonical lookup table from numeric tips back to standardized labels.

### `metadata.csv`

Per-tree summary table with fields including:

- standardized filename
- group
- study
- tip count
- dated flag

This file is the quickest way to audit which representative trees are present in the archive.

### `full_mapping.csv`

Expanded mapping table linking:

- numeric ID
- original source label
- standardized label
- source group

This is the traceability bridge between raw labels and the public numeric-tip archive.

---

## Important caveats

This folder now contains only archive components. Legacy figure byproducts (`figure_coverage.png`, `figure_cumulative.png`, `figure_dated.png`) have been moved to `phylo_trees/data/`.

**Fish tip-count note:** The standardized `fish.nwk` contains the 11,638-tip dated backbone from Rabosky 2018, which is the tree used for crown-age and LTT analyses. The manuscript coverage figures cite 31,516 tips, which refers to the FishTreeOfLife full chronogram with taxonomic imputation. The provenance entry documents both uses.

---

## Relationship to other layers

- Raw/source trees live under `phylo_trees/raw/`
- Working classification subsets live in `phylo_trees/curated/dated/` and `phylo_trees/curated/undated/`
- Derived analysis products live in `phylo_trees/data/`
- Manuscript-facing provenance lives in `phylo-datasets-paper/data_provenance.csv`

The archive should be auditable back to the raw/source layer and forward to the manuscript layer.

---

## Audit checklist for this folder

When auditing the standardized archive, verify that:

1. every intended representative dataset has the expected `.nwk` file in `trees/`
2. `metadata.csv` matches the representative tree actually used in the paper
3. dated flags are consistent with manuscript claims and derived figure inputs
4. `dictionary.csv` and `full_mapping.csv` have plausible sizes and no obvious label-cleaning failures
5. non-archive byproducts are either removed or explicitly classified as legacy outputs
