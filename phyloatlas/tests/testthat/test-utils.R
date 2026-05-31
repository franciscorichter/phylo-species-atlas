test_that(".provenance_key() identity-maps groups that already exist in provenance", {
  res <- phyloatlas:::.provenance_key(
    md_groups   = c("mammals", "fungi"),
    prov_groups = c("mammals", "fungi", "birds_mctavish")
  )
  expect_identical(res, c("mammals", "fungi"))
})

test_that(".provenance_key() rewrites 'birds' to 'birds_mctavish'", {
  res <- phyloatlas:::.provenance_key(
    md_groups   = "birds",
    prov_groups = c("birds_mctavish")
  )
  expect_identical(res, "birds_mctavish")
})

test_that(".provenance_key() collapses any condamine_* prefix to 'condamine'", {
  res <- phyloatlas:::.provenance_key(
    md_groups   = c("condamine_Vangidae", "condamine_Foo", "condamine_"),
    prov_groups = "condamine"
  )
  expect_identical(res, c("condamine", "condamine", "condamine"))
})

test_that(".provenance_key() returns NA_character_ for groups absent from provenance", {
  res <- phyloatlas:::.provenance_key(
    md_groups   = c("mammals", "unknown_group"),
    prov_groups = "mammals"
  )
  expect_identical(res, c("mammals", NA_character_))
})

test_that(".provenance_key() preserves input length and returns a character vector", {
  res <- phyloatlas:::.provenance_key(
    md_groups   = character(0),
    prov_groups = "mammals"
  )
  expect_type(res, "character")
  expect_length(res, 0L)

  res2 <- phyloatlas:::.provenance_key(
    md_groups   = c("mammals", "birds", "condamine_Vangidae", "orphan"),
    prov_groups = c("mammals", "birds_mctavish", "condamine")
  )
  expect_length(res2, 4L)
  expect_identical(res2, c("mammals", "birds_mctavish", "condamine", NA_character_))
})
