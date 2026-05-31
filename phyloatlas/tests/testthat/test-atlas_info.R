test_that("atlas_info() returns a one-row data frame for a known tree", {
  info <- atlas_info("mammals")
  expect_s3_class(info, "data.frame")
  expect_equal(nrow(info), 1L)
  expect_equal(info$name, "mammals")
  expect_equal(info$journal, "Nature")
  expect_identical(rownames(info), "1")
})

test_that("atlas_info() returns NULL with a warning for an unknown name", {
  expect_warning(out <- atlas_info("does_not_exist"), "No tree named")
  expect_null(out)
})

test_that("atlas_info() strips a trailing .nwk from the requested name", {
  expect_equal(atlas_info("birds.nwk")$name, "birds")
})

test_that("atlas_info() validates its input", {
  expect_error(atlas_info(character(0)))
  expect_error(atlas_info(c("a", "b")))
  expect_error(atlas_info(123))
})
