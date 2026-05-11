# phyloatlas

R access to the [Phylo-Species Atlas](https://github.com/franciscorichter/phylo-species-atlas) — 49 datasets, 641,763 species, standardized Newick format.

## Install

```r
# install.packages("pak")
pak::pkg_install("franciscorichter/phylo-species-atlas/phyloatlas")
```

## Use

```r
library(phyloatlas)

# Browse available trees
trees <- list_trees()
head(trees)

# Load one
tree <- load_atlas_tree("mammals")
plot(tree, show.tip.label = FALSE)

# Provenance for a single tree
atlas_info("birds")
```

Trees larger than ~10k tips download faster with `resolve_labels = FALSE`, which skips the 18 MB species dictionary:

```r
tree <- load_atlas_tree("seed_plants", resolve_labels = FALSE)
```

## Functions

| Function | Purpose |
|----------|---------|
| `load_atlas_tree(name)` | Load a tree by name, with standardized labels |
| `list_trees()` | Data frame of every tree + provenance |
| `atlas_info(name)` | One-row metadata for a single tree |
| `atlas_clear_cache()` | Force re-download of cached metadata/dictionary |

## Notes

Dictionary and metadata are cached in memory once per R session. To point at a fork or local mirror:

```r
options(phyloatlas.base_url = "https://raw.githubusercontent.com/yourfork/phylo-species-atlas/main")
atlas_clear_cache()
```
