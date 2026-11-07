# cran-comments

## Update

This is an update of phyloatlas 0.1.0 (on CRAN since 2026-06-04). One change:
`load_atlas_tree()` stops with an error that gives the address of the source
when the requested tree is one the atlas no longer redistributes, because its
source states no licence (data release v1.1.3 of the atlas). The list of such
trees is read from the atlas's `registers/withheld_trees.csv`, or from a copy
bundled in `inst/extdata/` when the atlas has none.

## Test environments

* local macOS aarch64, R 4.4.0: `R CMD check --as-cran` on the built tarball.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE: unable to verify current time.
  The check machine could not reach a time server; this is not a property of
  the package.

## Network handling

As in 0.1.0: the package fetches data over HTTPS from the atlas's public
GitHub repository. Network-dependent examples are in `\donttest{}` and guarded
with `try(..., silent = TRUE)`; the tests are offline (a `file://` mirror under
`tests/testthat/fixtures/`).

## Reverse dependencies

None.
