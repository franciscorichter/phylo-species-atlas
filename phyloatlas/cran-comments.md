# cran-comments

## Test environments

* local macOS aarch64, R 4.4.0
* (pre-submission) win-builder (devel, release)
* (pre-submission) R-hub multi-platform check

## R CMD check results

0 errors | 0 warnings | 1 expected NOTE on first submission:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Francisco Richter <richtf@usi.ch>'
New submission
```

This is the standard "New submission" NOTE expected for any first-time
CRAN submission. No other NOTEs or WARNINGs.

## Package purpose

`phyloatlas` is a thin R client for the Phylo-Species Atlas
(<https://github.com/franciscorichter/phylo-species-atlas>), a curated
collection of standardized empirical species-level phylogenetic trees.
All four exported functions fetch data over HTTPS from the atlas's
public GitHub repository.

## Network handling

The package fetches data from <https://raw.githubusercontent.com/franciscorichter/phylo-species-atlas/>.
Network-dependent examples are wrapped in `\donttest{}` (per the
Writing R Extensions recommendation over `\dontrun{}`) and additionally
guarded with `try(..., silent = TRUE)` so they do not raise errors
if the network is unavailable. Examples for `load_atlas_tree()` also
include an offline path that loads a small bundled tree from
`inst/extdata/`, so the function has at least one runnable example
under `R CMD check`.

Tests are entirely offline: the test suite sets
`options(phyloatlas.base_url = ...)` to a `file://` URL pointing at a
miniature mirror of the atlas under `tests/testthat/fixtures/`, so no
network access is required during `R CMD check`.

## Downstream dependencies

None — this is a new package.
