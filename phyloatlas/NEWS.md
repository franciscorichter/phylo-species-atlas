# phyloatlas 0.1.1

* `load_atlas_tree()` stops with an error that gives the address of the source
  when the requested tree is one the atlas does not redistribute because its
  source states no licence (data release v1.1.3: cephalopods, frogs,
  gastropods, insects, primates, seed_plants, sharks, sponges). The list is
  read from the atlas's `registers/withheld_trees.csv`; when the atlas named by
  `phyloatlas.base_url` has no such register, the copy bundled with the
  package (`inst/extdata/withheld_trees.csv`) is used.

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
