test_that("real atlas base URL still serves metadata.csv", {
  testthat::skip_on_cran()
  testthat::skip_if_offline("raw.githubusercontent.com")

  withr::with_options(list(phyloatlas.base_url = NULL), {
    atlas_clear_cache()
    out <- list_trees()
    expect_s3_class(out, "data.frame")
    expect_true(nrow(out) > 50)  # atlas has 76 trees as of 2026-05
    expect_true("name" %in% names(out))
  })
  atlas_clear_cache()
})
