# phyloatlas 0.1.0

* Initial CRAN release.
* Four user-facing functions:
  * `load_atlas_tree()` — fetch a standardized Newick tree from the Phylo-Species Atlas, optionally resolving integer tip IDs to species names via the shared dictionary.
  * `list_trees()` — data frame of every tree in the atlas with provenance fields (study, year, journal, DOI, coverage).
  * `atlas_info()` — single-row provenance for a tree by name.
  * `atlas_clear_cache()` — drop the in-session dictionary/metadata cache.
* In-memory caching of the dictionary and metadata files so repeated calls in a session don't re-download.
* Fork / mirror support via the `phyloatlas.base_url` option.
* Bundled offline demo data under `inst/extdata/` so examples and tests run without network access.
