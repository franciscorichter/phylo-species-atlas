# atlas_coverage -- the whole coverage layer of the phylogenetic atlas.
#
# Every published coverage figure in the paper is ONE call of this function
# against Data/coverage_basis.csv.  Numerator and denominator are column sums
# of the SAME 62-row file, which is what makes a phantom denominator
# impossible.
#
#   basis <- load_basis("Data/coverage_basis.csv")
#   euk   <- function(d) d$is_prokaryotic == "no"
#   atlas_coverage(basis, euk)                                    # C0
#   atlas_coverage(basis, euk, units_col = "c1_units")            # C1
#   atlas_coverage(basis, function(d) euk(d) & d$is_chronogram == "yes",
#                  units_col = "c1_units", denominator_predicate = euk)  # C1d
#
# c1_units is the number of distinct accepted names among the molecular tips
# (c1_basis = enumerated), c0_units where every tip is molecular (all_tips),
# and the expectation c0_units * f_mol only where neither holds.

load_basis <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

atlas_coverage <- function(basis, predicate = NULL, weight_col = NULL,
                           denominator_predicate = NULL,
                           units_col = "c0_units") {
  if (is.null(predicate)) predicate <- function(d) rep(TRUE, nrow(d))
  if (is.null(denominator_predicate)) denominator_predicate <- predicate
  keep <- predicate(basis)
  keep[is.na(keep)] <- FALSE
  w <- if (is.null(weight_col)) rep(1, nrow(basis)) else
    suppressWarnings(as.numeric(basis[[weight_col]]))
  w[is.na(w)] <- 0
  tips <- suppressWarnings(as.numeric(basis[[units_col]]))
  tips[is.na(tips)] <- 0
  num <- sum(tips[keep] * w[keep])
  dkeep <- denominator_predicate(basis)
  dkeep[is.na(dkeep)] <- FALSE
  desc <- suppressWarnings(as.numeric(basis$described))
  desc[is.na(desc)] <- 0
  den <- sum(desc[dkeep])
  list(numerator = num, denominator = den, n_partitions = sum(keep),
       percent = if (den > 0) 100 * num / den else NA_real_)
}
