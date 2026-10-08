# taxon_addition_dates.R -- does a taxon-addition tree set vary the dates of
# its sequenced species?  (decision D55, level C1u)
#
# Each tree of the fish all-taxon-assembled set (Rabosky et al. 2018,
# actinopt_full.trees) is pruned to the tips of the source's dated tree of
# sequenced species (actinopt_12k_treePL.tre), and its sorted branching times
# are compared with those of that tree.  Output: one row per set in
# Data/taxon_addition_set_check.csv.  The raw files are not redistributed
# (fishtreeoflife.org).
#
#   Rscript Paper/code/analysis/taxon_addition_dates.R   (from the repository root)
suppressMessages(library(ape))
d <- file.path("Data", "raw", "fish")
mol <- read.tree(file.path(d, "actinopt_12k_treePL.tre"))
bt0 <- sort(branching.times(mol))
trs <- read.tree(file.path(d, "actinopt_full.trees"))
res <- vapply(seq_along(trs), function(i) {
  t <- trs[[i]]
  p <- keep.tip(t, intersect(t$tip.label, mol$tip.label))
  bt <- sort(branching.times(p))
  c(Ntip(p), if (length(bt) == length(bt0)) max(abs(bt - bt0)) else NA_real_)
}, numeric(2))
out <- data.frame(
  partition = "Fish",
  set = "raw/fish/actinopt_full.trees",
  reference = "raw/fish/actinopt_12k_treePL.tre",
  n_trees = length(trs),
  n_sequenced_tips = Ntip(mol),
  n_tips_kept_min = min(res[1, ]),
  n_tips_kept_max = max(res[1, ]),
  max_abs_branching_time_diff_myr = sprintf("%.6g", max(res[2, ])))
write.csv(out, file.path("Data", "taxon_addition_set_check.csv"),
          row.names = FALSE, quote = FALSE)
print(out)
