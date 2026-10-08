# shark_set_comparison.R -- why the two shark tree sets give different fits
# (round 1 of 2026-10-09, finding F4-23)
#
# Stein et al. 2018 archive, for the 10-fossil scenario, (a) 500 trees of the
# 610 sequenced species, each the fixed RAxML topology dated in treePL under
# one root age drawn from the calibration prior, and (b) 10,000 taxon-complete
# trees: 582 species added to each of the 500 dated trees and the result
# resolved with Polytomy Resolver (BEAST, prior only), 20 draws per dated tree.
# Pruned to the 610 sequenced tips, the restriction of the worked example binds
# in far fewer trees of (b) than of (a).  This script measures whether the
# taxon addition and resolution steps change the topology and the node ages of
# the 610 sequenced species.
#
# Trees compared: the 100 trees of each set that worked_example_sets.R fits
# (indices round(seq(1, N, length.out = 100)), file order).  Each tree of (b),
# pruned to the 610 tips, is compared with every tree of (a): its clades are
# matched by tip set, and the tree of (a) whose clade ages differ least (median
# relative difference over shared clades) is its best match.
#
# Output: Data/shark_set_comparison.csv (one row).  The raw files are not
# redistributed (data.vertlife.org/sharktree).
#
#   Rscript Paper/code/analysis/shark_set_comparison.R   (from the repository root)
suppressMessages(library(ape))
f_mol <- file.path("Data", "raw", "sharks", "Chond.610sp.10Cal.500TreeSet.tre")
f_cmp <- file.path("Data", "raw", "revision-2026-10", "sharks", "ten_calibration",
                   "10.cal.tree.nex")
REL_TOL <- 0.01      # clade ages within 1% count as equal
YOUNG_MYR <- 5       # nodes younger than this count as young
N_SELECT <- 100

mol_all <- read.tree(f_mol)
labs <- sort(mol_all[[1]]$tip.label)
sel <- function(n) unique(round(seq(1, n, length.out = N_SELECT)))

# (b): read only the selected trees of the 357 MB nexus file
ln <- readLines(f_cmp)
is_tree <- grepl("^\\s*TREE ", ln)
n_cmp <- sum(is_tree)
keep <- which(is_tree)[sel(n_cmp)]
hdr <- ln[seq_len(which(is_tree)[1] - 1)]
tmp <- tempfile(fileext = ".nex")
writeLines(c(hdr, ln[keep], "END;"), tmp)
rm(ln)
cmp <- read.nexus(tmp)
unlink(tmp)

clades <- function(t) {
  t <- keep.tip(t, labs)
  pp <- prop.part(t)
  l <- attr(pp, "labels")
  k <- vapply(pp, function(x) paste(sort(l[x]), collapse = "|"), "")
  setNames(branching.times(t)[as.character(Ntip(t) + seq_along(pp))], k)
}
topo <- function(t) { t$edge.length <- NULL; write.tree(t) }
n_topo_mol <- length(unique(vapply(mol_all, topo, "")))
M <- vapply(mol_all, clades, numeric(Ntip(mol_all[[1]]) - 1))  # fixed topology
kM <- rownames(M)
gen <- function(x) sub("_.*", "", x)

res <- t(vapply(cmp, function(t) {
  h <- clades(t)
  sh <- intersect(names(h), kM)
  d <- apply(M[sh, , drop = FALSE], 2, function(m) median(abs(m - h[sh]) / m))
  j <- which.min(d)
  rel <- abs(M[sh, j] - h[sh]) / M[sh, j]
  absent <- setdiff(kM, names(h))
  within <- vapply(strsplit(absent, "|", fixed = TRUE),
                   function(x) length(unique(gen(x))) == 1, logical(1))
  c(shared = length(sh), absent = length(absent), absent_within_genus = sum(within),
    equal_age = sum(rel <= REL_TOL), young = sum(h < YOUNG_MYR),
    gamma = gammaStat(keep.tip(t, labs)),
    root_diff = abs(max(h) - max(M[, j])))
}, numeric(7)))
mol <- mol_all[sel(length(mol_all))]
young_mol <- vapply(mol, function(t) sum(branching.times(t) < YOUNG_MYR), 0)
gamma_mol <- vapply(mol, gammaStat, 0)

f3 <- function(x) sprintf("%.3f", x)
f2 <- function(x) sprintf("%.2f", x)
out <- data.frame(
  partition = "Sharks",
  complete_set = sub("^Data/", "", f_cmp),
  sequenced_set = sub("^Data/", "", f_mol),
  n_complete_trees = n_cmp, n_sequenced_trees = length(mol_all),
  n_compared = nrow(res), n_internal_nodes = length(kM),
  sequenced_set_topologies = n_topo_mol,
  pruned_topology_equal = sum(res[, "absent"] == 0),
  absent_clades_median = median(res[, "absent"]),
  absent_clades_min = min(res[, "absent"]), absent_clades_max = max(res[, "absent"]),
  absent_within_genus_median = median(res[, "absent_within_genus"]),
  equal_age_tol = REL_TOL,
  equal_age_median = median(res[, "equal_age"]),
  equal_age_min = min(res[, "equal_age"]), equal_age_max = max(res[, "equal_age"]),
  root_diff_median_myr = f3(median(res[, "root_diff"])),
  young_myr = YOUNG_MYR,
  young_complete_median = median(res[, "young"]),
  young_sequenced_median = median(young_mol),
  gamma_complete_median = f2(median(res[, "gamma"])),
  gamma_sequenced_median = f2(median(gamma_mol)))
write.csv(out, file.path("Data", "shark_set_comparison.csv"),
          row.names = FALSE, quote = FALSE)
print(t(out))
