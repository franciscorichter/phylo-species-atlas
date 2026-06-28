# Phylogenetic Species Atlas

A curated, standardized collection of empirical species-level phylogenies across the tree of life, with per-tree provenance and an archival-uncertainty audit.

## Overview

This atlas surveys publicly available, downloadable species-level empirical phylogenies across all major branches of life. It accompanies the paper:

> **Where the tree of life is empirically resolved, and where it is not: an open atlas of species-level phylogenies and their archival uncertainty**
>
> Francisco Richter, Università della Svizzera italiana (USI)
> Zenodo (concept DOI): https://doi.org/10.5281/zenodo.20127157

### At a glance

| Item | Value |
|---|---|
| Standardized trees | **264** |
| Atlas partitions | **62** (26 with shipped canonical + 36 unrepresented) |
| Source datasets | **49** |
| Taxonomic groups | **47** |
| Decomposition | 218 Condamine family-level + 25 further partition-canonical + 21 sub-clade/reference |
| Global dictionary | **637,619** standardized labels (493,978 eukaryotic + 143,641 prokaryotic) |
| Time-calibrated | **248 of 264 (94%)** |
| Eukaryotic coverage (headline C0) | **23.7%** of ~2.1 M described |
| Sensitivity bound C1 (strict molecular-only on TACT-imputed trees) | **~10%** |
| Sensitivity bound C2 (strictest non-imputed, non-OTL-pipeline) | **~5.4%** |
| Archival uncertainty preserved | **7 of 24** non-Condamine dated source trees |
| Catalogue of Life 2025 reconciliation | included |
| Licence | CC-BY-4.0 (see `LICENSE`) |
| Version | v1.0.5 (2026-06-13) |

### Three-category tree classification

Trees in the atlas are classified explicitly into three categories:

1. **Direct empirical phylogenies** — inferred from molecular or genomic data in a primary study (e.g., Upham et al. 2019 mammals; Stein et al. 2018 sharks).
2. **Peer-reviewed empirical syntheses** — clade-bounded trees that combine published input phylogenies, sometimes augmented with taxonomic placement, released as standalone peer-reviewed deliverables (e.g., McTavish et al. 2025 birds; Smith & Brown 2018 seed plants; Chesters 2017 insects). May use Open Tree of Life pipeline software, but the resulting tree is published, citable, and clade-bounded.
3. **Excluded**: the live ~2.4M-tip Open Tree of Life synthesis product (its tips are placed predominantly by taxonomic descent and it lacks a clade-bounded primary methodology paper).

Per-tree provenance — molecular-tip fraction, Open Tree pipeline involvement, dating regime, archival-uncertainty status — is documented in Supplementary Table S5 of the manuscript and in the per-partition `info.yaml` files at `Data/site/data/partitions/<slug>/`.

## Repository structure

```
.
├── README.md
├── standardized/
│   ├── trees/              264 Newick files (numeric tip labels)
│   ├── dictionary.csv      Standardized-label dictionary: ID → standardized name (637,619 entries)
│   ├── full_mapping.csv    Full traceability: ID, original label, standardized name, group
│   ├── col_mapping.csv     Catalogue of Life 2025 name mapping (637,619 entries)
│   └── metadata.csv        Per-tree metadata: filename, group, study, tips, dated
└── data_provenance.csv     Complete metadata for all 49 datasets
```

## Standardized tree files

All trees are distributed in Newick format with **numeric tip labels** (integers 1--637,619). To recover species names, join the numeric labels against `dictionary.csv`.

### Load any tree from R: the `phyloatlas` package

```r
# install.packages("pak")
pak::pkg_install("franciscorichter/phylo-species-atlas/phyloatlas")

library(phyloatlas)

list_trees()                                 # all 264 trees + provenance
atlas_info("birds")                          # metadata for one tree
tree <- load_atlas_tree("mammals")           # Upham 2019 (5,912 species)
tree <- load_atlas_tree("birds")             # McTavish 2025 (10,824 species)
tree <- load_atlas_tree("seed_plants")       # Smith & Brown 2018 (356,305 species)
tree <- load_atlas_tree("condamine_Vangidae")# Condamine family tree (21 species)
plot(tree)
```

For very large trees, skip the species-name lookup (keeps integer IDs, avoids the 18 MB dictionary download):

```r
tree <- load_atlas_tree("seed_plants", resolve_labels = FALSE)
```

See [`phyloatlas/README.md`](phyloatlas/README.md) for the full API. If you'd rather not install a package, the equivalent ~5-line snippet using just `ape` works too:

```r
library(ape)
base <- "https://raw.githubusercontent.com/franciscorichter/phylo-species-atlas/main"
tree <- read.tree(paste0(base, "/standardized/trees/mammals.nwk"))
dict <- read.csv(paste0(base, "/standardized/dictionary.csv"))
tree$tip.label <- dict$standardized_name[match(as.integer(tree$tip.label), dict$id)]
```

### Load from a local clone

```r
library(ape)

tree <- read.tree("standardized/trees/mammals.nwk")
dict <- read.csv("standardized/dictionary.csv")
tree$tip.label <- dict$standardized_name[match(as.integer(tree$tip.label), dict$id)]
```

## Species dictionary

`dictionary.csv` contains two columns:

| Column | Description |
|--------|-------------|
| `id` | Sequential integer (1--637,619) |
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

Match rates: **83.2% of matchable eukaryotic species** mapped to CoL (352K exact + 58K via synonym resolution). The 143K GTDB prokaryote entries use genome accession IDs not present in CoL.

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

## License

The standardized tree files, dictionary, and metadata are released for academic use. Individual source trees retain the licenses of their original publications. See `data_provenance.csv` for DOIs and original data sources.
