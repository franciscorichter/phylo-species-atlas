test_that("atlas_clear_cache() empties the internal cache", {
  atlas_clear_cache()
  expect_length(ls(phyloatlas:::.cache), 0L)
  invisible(list_trees())
  expect_true(length(ls(phyloatlas:::.cache)) > 0L)
  expect_true(all(c("metadata", "provenance") %in% ls(phyloatlas:::.cache)))
  res <- atlas_clear_cache()
  expect_null(res)
  expect_length(ls(phyloatlas:::.cache), 0L)
})

test_that("After clear, list_trees() re-reads from the (possibly new) base URL", {
  # Build a second fixture root that has ONLY mammals.
  tmp <- file.path(tempdir(), "alt-atlas-root")
  dir.create(file.path(tmp, "standardized", "trees"), recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c("filename,group,study,ntips,dated",
      "mammals.nwk,mammals,Upham et al. 2019,4,TRUE"),
    file.path(tmp, "standardized", "metadata.csv")
  )
  file.copy(
    testthat::test_path("fixtures", "atlas-root", "standardized", "dictionary.csv"),
    file.path(tmp, "standardized", "dictionary.csv"),
    overwrite = TRUE
  )
  file.copy(
    testthat::test_path("fixtures", "atlas-root", "data_provenance.csv"),
    file.path(tmp, "data_provenance.csv"),
    overwrite = TRUE
  )
  file.copy(
    testthat::test_path("fixtures", "atlas-root", "standardized", "trees", "mammals.nwk"),
    file.path(tmp, "standardized", "trees", "mammals.nwk"),
    overwrite = TRUE
  )

  withr::with_options(
    list(phyloatlas.base_url = paste0("file://", utils::URLencode(normalizePath(tmp, winslash = "/")))),
    {
      atlas_clear_cache()  # critical - without this we'd see the stale cache
      out <- list_trees()
      expect_equal(nrow(out), 1L)
      expect_equal(out$name, "mammals")
    }
  )
})
