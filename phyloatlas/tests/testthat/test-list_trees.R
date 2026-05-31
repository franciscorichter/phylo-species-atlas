test_that("list_trees() returns a data frame ordered by group with one row per tree", {
  out <- list_trees()
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  # Sorted by group (birds, condamine_Vangidae, mammals alphabetically)
  expect_identical(out$group, c("birds", "condamine_Vangidae", "mammals"))
  expect_identical(rownames(out), as.character(seq_len(nrow(out))))
})

test_that("list_trees() exposes the documented columns and only those", {
  out <- list_trees()
  expected <- c(
    "name", "group", "study", "ntips", "dated",
    "year", "journal", "doi", "crown_ma",
    "described_species", "coverage_pct",
    "data_source", "download_url", "methods_brief", "notes",
    "study_full"
  )
  expect_setequal(names(out), expected)
  # column ordering matches the explicit cols vector in list_trees()
  expect_identical(names(out), expected)
})

test_that("list_trees() joins provenance correctly across all three key-mapping branches", {
  out <- list_trees()
  # mammals -> direct
  mammals <- out[out$name == "mammals", ]
  expect_equal(mammals$year, 2019)
  expect_equal(mammals$journal, "Nature")
  # birds -> renamed to birds_mctavish in provenance
  birds <- out[out$name == "birds", ]
  expect_equal(birds$year, 2025)
  # condamine_Vangidae -> collapsed to condamine in provenance
  vang <- out[out$name == "condamine_Vangidae", ]
  expect_equal(vang$year, 2019)
  expect_equal(vang$data_source, "Dryad")
})

test_that("list_trees() renames provenance 'study' to 'study_full' to avoid clashing with metadata 'study'", {
  out <- list_trees()
  expect_true("study" %in% names(out))
  expect_true("study_full" %in% names(out))
  expect_equal(out$study[out$name == "mammals"], "Upham et al. 2019")
  expect_equal(out$study_full[out$name == "mammals"], "Upham 2019 long form")
})
