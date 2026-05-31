# Compute an absolute, platform-safe file:// URL for the bundled atlas mirror.
.fixture_root <- function() {
  root <- testthat::test_path("fixtures", "atlas-root")
  paste0("file://", utils::URLencode(normalizePath(root, winslash = "/", mustWork = TRUE)))
}

# Activate the fixture base URL for the whole test session.
old_opts <- options(phyloatlas.base_url = .fixture_root())

# Wipe the package cache before every test so caching tests are deterministic
# and option changes (e.g. test-atlas_clear_cache.R) take effect immediately.
testthat::setup({
  if (exists("atlas_clear_cache", mode = "function")) atlas_clear_cache()
})

testthat::teardown({
  if (exists("atlas_clear_cache", mode = "function")) atlas_clear_cache()
})
