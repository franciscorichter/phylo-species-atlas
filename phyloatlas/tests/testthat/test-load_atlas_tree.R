test_that("load_atlas_tree() returns a phylo object and resolves dictionary ids", {
  tr <- load_atlas_tree("mammals")
  expect_s3_class(tr, "phylo")
  expect_equal(length(tr$tip.label), 4L)
  expect_true(all(c("Homo_sapiens", "Mus_musculus", "Pan_troglodytes") %in% tr$tip.label))
})

test_that("load_atlas_tree() leaves tip labels untouched when the integer is not in the dictionary", {
  tr <- load_atlas_tree("mammals")
  # id 999 has no dictionary row -> stays as the original string "999"
  expect_true("999" %in% tr$tip.label)
})

test_that("load_atlas_tree(resolve_labels = FALSE) returns raw integer-string tip labels and does not touch the dictionary cache", {
  atlas_clear_cache()
  tr <- load_atlas_tree("birds", resolve_labels = FALSE)
  expect_setequal(tr$tip.label, c("4", "5"))
  # dictionary must NOT have been loaded by this code path
  expect_false("dictionary" %in% ls(phyloatlas:::.cache))
})

test_that("load_atlas_tree() strips a trailing .nwk from the requested name", {
  tr1 <- load_atlas_tree("birds")
  tr2 <- load_atlas_tree("birds.nwk")
  expect_identical(sort(tr1$tip.label), sort(tr2$tip.label))
})

test_that("load_atlas_tree() errors with a helpful message on an unknown tree name", {
  expect_error(
    load_atlas_tree("does_not_exist"),
    regexp = "Could not load tree 'does_not_exist'"
  )
})

test_that("load_atlas_tree() validates its input", {
  expect_error(load_atlas_tree(character(0)))
  expect_error(load_atlas_tree(c("a", "b")))
  expect_error(load_atlas_tree(""))
  expect_error(load_atlas_tree(42))
})
