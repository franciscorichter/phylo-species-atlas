#!/usr/bin/env Rscript
# worked_example_sets.R -- the worked example repeated over the tree sets the
# sources archive (round F2: findings F2-04, F2-05, F2-35).
#
# QUESTION.  worked_example_bd.R fits one summary tree per group.  Three of
# those trees are consensus summaries whose polytomies the likelihood reads as
# simultaneous splits, and the bird summary tree is one of many date
# assignments.  Where the source archives the trees it summarised, the fits
# are repeated on a fixed, evenly spaced sample of them, so the summary-tree
# value can be read as one point of a distribution.
#
# SETS (Data/worked_example_tree_sets.csv, one row per set):
#   tips = molecular   each tree is pruned to the molecular tips of the
#                      reference tree (Data/molecular_tips.csv), matched by
#                      their source labels (Data/standardized/full_mapping.csv,
#                      original label or standardized name); the full tree
#                      (outgroups removed) is fitted too.
#   tips = all         the source's own tree of its molecular species: every
#                      tip (outgroups removed) is molecular.
#   Selection: the trees of the set in file order (files sorted by name),
#   indices round(seq(1, N, length.out = n_select)).  No tree is rescaled:
#   each carries its own dates.
#
# FITS per tree, with the likelihood, optimiser and conditioning of
# worked_example_bd.R (copied below, unchanged):
#   matching   pruned tree, rho_m = its tips / described species
#   full rho   pruned tree, rho_f = reference-tree tips / described species
#   full tree  (tips = molecular only) the full tree, rho = its tips /
#              described species
#   LR = 2 (logL_matching - logL_full_rho); the restriction is rejected at the
#   5% level when LR > qchisq(0.90, 1) (Self & Liang 1987), as in the main fits.
#
# OUTPUT  Paper/inputs/worked_example_bd_sets.csv   one row per tree
#
# USAGE   Rscript Paper/code/analysis/worked_example_sets.R [REPO_ROOT]

suppressPackageStartupMessages({ library(ape) })

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(dirname(f), "..", "..", ".."))
}
P <- function(...) file.path(root, ...)

EPS_MAX <- 1 - 1e-9
ZERO_BRANCH_TOL <- 1e-6

basis <- read.csv(P("Data", "coverage_basis.csv"), stringsAsFactors = FALSE)
moltips <- read.csv(P("Data", "molecular_tips.csv"), stringsAsFactors = FALSE,
                    colClasses = "character")
fmap <- read.csv(P("Data", "standardized", "full_mapping.csv"),
                 stringsAsFactors = FALSE, colClasses = "character")
sets <- read.csv(P("Data", "worked_example_tree_sets.csv"),
                 stringsAsFactors = FALSE)
# WE_SETS_ONLY="Birds,Sharks" restricts a run to some sets and writes to
# WE_SETS_OUT (a partial file).  Where the sets are held on different hosts,
# the committed file is the concatenation, in register order, of partial runs
# of this script; each row depends only on its own tree.
only <- Sys.getenv("WE_SETS_ONLY")
OUT <- P("Paper", "inputs", "worked_example_bd_sets.csv")
if (nzchar(only)) {
  sets <- sets[sets$partition %in% strsplit(only, ",")[[1]], ]
  OUT <- Sys.getenv("WE_SETS_OUT",
                    file.path(tempdir(), "worked_example_bd_sets_partial.csv"))
}

# ---- likelihood (identical to worked_example_bd.R) -------------------------
branching_events <- function(tr) {
  n <- Ntip(tr)
  depth <- node.depth.edgelength(tr)
  h <- max(depth[1:n])
  ages <- h - depth[(n + 1):(n + tr$Nnode)]
  kids <- tabulate(tr$edge[, 1], nbins = n + tr$Nnode)[(n + 1):(n + tr$Nnode)]
  ev <- rep(ages, kids - 1)
  sort(pmax(ev, 0), decreasing = TRUE)
}
polytomy_events <- function(tr) {
  t <- di2multi(tr, tol = ZERO_BRANCH_TOL)
  n <- Ntip(t)
  kids <- tabulate(t$edge[, 1], nbins = n + t$Nnode)[(n + 1):(n + t$Nnode)]
  sum((kids - 1L)[kids > 2L])
}
loglik <- function(r, eps, x, rho) {
  lam <- r / (1 - eps); mu <- eps * lam
  n <- length(x) + 1
  e <- exp(-r * x)
  den <- rho * lam + (lam * (1 - rho) - mu) * e
  lp1 <- log(rho) + 2 * log(r) - r * x - 2 * log(den)
  l1p0 <- log(rho) + log(r) - log(den[1])
  (n - 2) * log(lam) + 2 * lp1[1] + sum(lp1[-1]) - 2 * l1p0
}
fit <- function(x, rho) {
  nll <- function(p) -loglik(exp(p[1]), plogis(p[2]) * EPS_MAX, x, rho)
  best <- NULL
  for (e0 in c(-4, 0, 2)) {
    r0 <- log((length(x) + 1) / x[1]) / x[1]
    o <- optim(c(log(max(r0, 1e-4)), e0), nll, method = "Nelder-Mead",
               control = list(maxit = 4000, reltol = 1e-12))
    if (is.null(best) || o$value < best$value) best <- o
  }
  list(r = exp(best$par[1]), eps = plogis(best$par[2]) * EPS_MAX,
       ll = -best$value)
}

# ---- reading the sets -------------------------------------------------------
tree_texts <- function(s) {
  files <- sort(Sys.glob(P("Data", s$source)))
  if (!length(files)) stop("no file matches ", s$source)
  if (s$format == "newick_files") {
    idx <- unique(round(seq(1, length(files), length.out = s$n_select)))
    return(list(idx = idx, text = vapply(files[idx], function(f)
      paste(readLines(f, warn = FALSE), collapse = ""), ""), translate = NULL))
  }
  if (s$format == "newick_lines") {
    lines <- unlist(lapply(files, function(f) {
      l <- readLines(f, warn = FALSE); l[grepl("^\\s*\\(", l)] }))
    idx <- unique(round(seq(1, length(lines), length.out = s$n_select)))
    return(list(idx = idx, text = lines[idx], translate = NULL))
  }
  if (s$format == "nexus_translate") {
    l <- readLines(files[1], warn = FALSE)
    a <- grep("^\\s*TRANSLATE\\s*$", l, ignore.case = TRUE)[1]
    b <- a + which(grepl(";\\s*$", l[(a + 1):length(l)]))[1]
    kv <- trimws(sub(",\\s*$", "", sub(";\\s*$", "", l[(a + 1):b])))
    kv <- kv[nzchar(kv)]
    tr <- setNames(sub("^\\S+\\s+", "", kv), sub("\\s.*$", "", kv))
    tl <- grep("^\\s*TREE\\s", l, ignore.case = TRUE, value = TRUE)
    tl <- sub("^[^(]*", "", tl)
    idx <- unique(round(seq(1, length(tl), length.out = s$n_select)))
    return(list(idx = idx, text = tl[idx], translate = tr))
  }
  stop("unknown format ", s$format)
}

rows <- list()
for (i in seq_len(nrow(sets))) {
  s <- sets[i, ]
  b <- basis[basis$partition == s$partition, ]
  group <- sub("\\.nwk$", "", b$canonical_tree)
  ids <- moltips$standardized_id[moltips$partition == s$partition]
  fm <- fmap[fmap$group == group & fmap$id %in% ids, ]
  # a source label can carry suffixes (Upham et al.: genus_epithet_FAMILY_ORDER)
  # that a set's labels omit, so the genus and epithet are matched as well
  mol_labels <- unique(c(fm$original_label, fm$standardized_name,
                         sub("^([^_]+_[^_]+).*$", "\\1", fm$original_label)))
  outside <- trimws(strsplit(s$outside, ";")[[1]])
  outside <- outside[nzchar(outside)]
  rho_f <- b$labels_unique / b$described
  tt <- tree_texts(s)
  cat(sprintf("%-10s %-14s %d trees\n", s$partition, s$set, length(tt$idx)))
  for (k in seq_along(tt$idx)) {
    tr <- read.tree(text = tt$text[k])
    if (!is.null(tt$translate))
      tr$tip.label <- unname(ifelse(tr$tip.label %in% names(tt$translate),
                                    tt$translate[tr$tip.label], tr$tip.label))
    tr <- drop.tip(tr, intersect(outside, tr$tip.label))
    full <- NULL
    if (s$tips == "molecular") {
      keep <- intersect(tr$tip.label, mol_labels)
      full <- tr
      tr <- keep.tip(tr, keep)
    }
    n <- Ntip(tr)
    rho_m <- n / b$described
    x <- branching_events(tr)
    mm <- fit(x, rho_m)
    mf <- fit(x, rho_f)
    out <- data.frame(
      partition = s$partition, set = s$set, tree_index = tt$idx[k],
      n_tips = n, root_height_ma = round(x[1], 4),
      n_events = length(x), n_polytomy_events = polytomy_events(tr),
      rho_matching = signif(rho_m, 8), rho_full = signif(rho_f, 8),
      r_matching = signif(mm$r, 6), epsilon_matching = signif(mm$eps, 6),
      loglik_matching = round(mm$ll, 4),
      r_full_rho = signif(mf$r, 6), epsilon_full_rho = signif(mf$eps, 6),
      loglik_full_rho = round(mf$ll, 4),
      lr_statistic = round(max(0, 2 * (mm$ll - mf$ll)), 4),
      n_tips_full_tree = NA_integer_, r_full_tree = NA_real_,
      epsilon_full_tree = NA_real_, loglik_full_tree = NA_real_,
      stringsAsFactors = FALSE)
    if (!is.null(full)) {
      xf <- branching_events(full)
      rho_t <- min(Ntip(full) / b$described, 1)
      mt <- fit(xf, rho_t)
      out$n_tips_full_tree <- Ntip(full)
      out$r_full_tree <- signif(mt$r, 6)
      out$epsilon_full_tree <- signif(mt$eps, 6)
      out$loglik_full_tree <- round(mt$ll, 4)
    }
    rows[[length(rows) + 1]] <- out
    if (k %% 20 == 0)
      cat(sprintf("  %d/%d  n=%d eps_m=%.3f LR=%.2f\n", k, length(tt$idx), n,
                  mm$eps, out$lr_statistic))
  }
}
res <- do.call(rbind, rows)
write.csv(res, OUT, row.names = FALSE, na = "")
cat("wrote", OUT, "(", nrow(res), "rows)\n")
