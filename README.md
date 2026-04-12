# Phylogenetic Species Atlas

A curated collection of empirical phylogenetic trees spanning the tree of life, with a standardized species dictionary linking 641,763 unique species across 264 tree files from 49 independent datasets.

## Overview

This atlas surveys the current state of publicly available, downloadable species-level molecular phylogenies across all major branches of life. It accompanies the paper:

> **The phylogenetic dark matter: a cross-kingdom atlas of species-level tree coverage**
>
> Francisco Richter

### Key numbers

- **49 datasets** from peer-reviewed studies (2010--2026)
- **641,763 unique species** in the standardized dictionary
- **264 Newick files** with numeric tip labels
- **~23% of described eukaryotic species** represented
- Coverage spans vertebrates, arthropods, plants, fungi, protists, bacteria, and archaea

## Repository structure

```
.
├── README.md
├── standardized/
│   ├── trees/              264 Newick files (numeric tip labels)
│   ├── dictionary.csv      Species dictionary: ID → standardized name (641,763 entries)
│   ├── full_mapping.csv    Full traceability: ID, original label, standardized name, group
│   ├── col_mapping.csv     Catalogue of Life 2025 name mapping (641,763 entries)
│   └── metadata.csv        Per-tree metadata: filename, group, study, tips, dated
├── data_provenance.csv     Complete metadata for all 49 datasets
└── code/
    └── figures/            R scripts to reproduce all manuscript figures
```

## Standardized tree files

All trees are distributed in Newick format with **numeric tip labels** (integers 1--641,763). To recover species names, join the numeric labels against `dictionary.csv`:

```r
library(ape)

# Read a tree
tree <- read.tree("standardized/trees/mammals.nwk")

# Load dictionary
dict <- read.csv("standardized/dictionary.csv")

# Map numeric tips back to species names
tree$tip.label <- dict$standardized_name[match(as.integer(tree$tip.label), dict$id)]
```

## Species dictionary

`dictionary.csv` contains two columns:

| Column | Description |
|--------|-------------|
| `id` | Sequential integer (1--641,763) |
| `standardized_name` | Standardized species label in `Genus_species` format |

Labels were standardized by stripping family/order annotations, voucher codes, and institutional suffixes. GTDB genome accessions (bacteria, archaea) are retained as-is.

## Catalogue of Life mapping

`col_mapping.csv` links each atlas entry to its Catalogue of Life 2025 accepted name:

| Column | Description |
|--------|-------------|
| `id` | Dictionary ID |
| `standardized_name` | Atlas species label |
| `col_id` | CoL accepted taxon ID (empty if no match) |
| `col_name` | CoL accepted scientific name |
| `col_kingdom` | Kingdom from CoL classification |
| `col_family` | Family from CoL classification |
| `match_type` | `exact`, `synonym`, `gtdb_only`, `genus_only`, `skip`, or `unmatched` |

Match rates: **83.6% of matchable eukaryotic species** mapped to CoL (352K exact + 58K via synonym resolution). The 143K GTDB prokaryote entries use genome accession IDs not present in CoL.

## Full mapping

`full_mapping.csv` provides complete traceability from original tree tip labels to standardized names:

| Column | Description |
|--------|-------------|
| `id` | Dictionary ID |
| `original_label` | Original tip label as it appeared in the source tree |
| `standardized_name` | Cleaned label |
| `group` | Taxonomic group (e.g., mammals, birds, insects) |

## Data provenance

`data_provenance.csv` provides complete metadata for each of the 49 datasets, including DOIs, download URLs, phylogenetic methods, and species count sources. See the Supplementary Information for formatted versions of this data (Tables S6 and S7).

## Taxonomic groups covered

**Vertebrates** (14): mammals, birds, amphibians, frogs, salamanders, fish, neotropical fish, squamates, sharks, turtles, primates, parrots, carnivora, cetaceans

**Arthropods** (8): insects (2 datasets), butterflies, bees, ants, hemiptera, termites, spiders, crustaceans

**Plants** (10): seed plants, ferns, conifers, bryophytes, grasses, orchids, asteraceae, solanaceae, cacti, eucalyptus

**Other animals** (8): corals, cephalopods, nematodes, bryozoa, bivalves, gastropods, sponges, echinoderms

**Microbes and protists** (4): diatoms, bacteria, archaea, fungi

**Cross-kingdom** (2): eukaryotes (backbone), TimeTree

**Condamine families** (218): tetrapod family-level dated trees from Condamine et al. 2019

## Figure reproduction

R scripts in `code/figures/` reproduce all manuscript figures. They read from `standardized/` and `data_provenance.csv`.

## Citation

If you use this atlas, please cite:

```
Richter, F. (2026). The phylogenetic dark matter: a cross-kingdom atlas of
species-level tree coverage. [preprint]
```

## License

The standardized tree files, dictionary, and metadata are released for academic use. Individual source trees retain the licenses of their original publications. See `data_provenance.csv` for DOIs and original data sources.
