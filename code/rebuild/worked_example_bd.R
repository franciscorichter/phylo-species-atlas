#!/usr/bin/env Rscript
# worked_example_bd.R -- what happens to a birth-death estimate when a coverage
# rung is used as the sampling fraction (finding S13, decisions D1, D23).
#
# QUESTION.  A user who fits a constant-rate birth-death model to a canonical
# chronogram must give the fraction rho of the clade's species the tree
# samples.  The coherent value is the tree's own tips over described species
# (D4 makes one tip one species).  A user who prunes the tree to its molecular
# tips (C1) but keeps the sampling fraction of the full tree (its tips over
# described species, the fraction a user finds quoted for the downloaded
# tree; decision D27) credits the pruned tree with species it does not
# contain: that is the error the ladder is meant to prevent.  The reconstructed process identifies r = lambda - mu almost
# independently of rho (Stadler 2009), so the error shows in lambda, mu and
# epsilon = mu / lambda, which are reported on the SAME tree under both rho.
#
# FITS, for every eligible partition (no selection; all are reported):
#   full       the canonical chronogram,              rho = n_tips / described
#   molecular  the same tree pruned to its molecular
#              tips (Data/molecular_tips.csv),         rho = n_tips / described
#   molecular  the same pruned tree,                  rho = full-tree tips /
#              described (rho_convention 'full_misapplied', D27)
# The full-versus-molecular comparison is a comparison of two different trees,
# not a rung effect: it shows what taxonomically placed tips do to a rate
# estimate (Rabosky 2015 Evolution 69:3207-3216; Chang, Rabosky & Alfaro 2020
# Syst. Biol. 69:602-611), and its direction is reported per partition.
#
# ELIGIBLE: eukaryotic canonical chronograms in absolute time whose molecular
# tips are enumerated and fewer than all tips (so the pruned tree differs).
#
# TIME SCALE (round F3, F3-04; decision D42).  Each tree is fitted with the
# branch lengths its source archived (time_scale_factor = 1).  Where the
# source summarises the crown over a tree set, the summary tree's root height
# (geometry_root_ma) differs from that crown age (crown_age_ma, D20); the
# ratio is reported, not applied.  Multiplying every branch length by c
# divides r, lambda and mu by c and adds -(n - 2) log c to every
# log-likelihood of the tree, so epsilon, the likelihood-ratio statistics and
# the boundary classification do not depend on it.
#
# MODEL.  Constant-rate birth-death with uniform random sampling of extant
# species with probability rho (Stadler 2009 J. Theor. Biol. 261:58-66, eq. 2;
# the Nee et al. 1994 likelihood with sampling), conditioned on the crown age
# and on both crown lineages leaving sampled descendants:
#     log L = (n-2) log(lambda) + 2 log p1(x1) + sum_{i>=2} log p1(x_i)
#             - 2 log(1 - p0(x1))
#     p1(t) = rho r^2 e^{-rt} / (rho lambda + (lambda(1-rho) - mu) e^{-rt})^2
#     1 - p0(t) = rho r / (rho lambda + (lambda(1-rho) - mu) e^{-rt})
# with r = lambda - mu and x1 > x2 >= ... the n-1 branching times.  A node
# with k children contributes k-1 branching events at its age (a polytomy is
# read as simultaneous splits).  Parameters are fitted as r > 0 and
# epsilon = mu / lambda in [0, 1); 95% intervals by profile likelihood
# (likelihood-ratio cutoff qchisq(0.95, 1) / 2, the same for every
# interval).  Where the MLE of epsilon is on its boundary (epsilon = 0) the
# interval is [0, upper]: 0 is never rejected (its statistic is 0), and every
# epsilon > 0 is an interior value whose statistic is chi2(1), so the cutoff
# that gives 95% coverage there is the ordinary one.  (A cutoff of
# qchisq(0.90, 1) / 2 applied only when the MLE is at 0 chooses the interval
# from the data and covers an interior epsilon with probability 0.90; it is
# not used.)  The column epsilon_ci says whether the MLE is at 0.
#
# THE IDENTITY.  The likelihood above depends on lambda and rho only through
# lambda * rho and r (Stadler 2009).  Changing rho to rho' with r fixed
# therefore leaves the likelihood unchanged when lambda' = lambda rho / rho'
# keeps mu' = lambda' - r >= 0: an interior fit at rho' is the fit at rho
# with lambda rescaled by rho / rho', and the log-likelihoods are equal.  The
# full-tree rho on the molecular tree is thus not evidence about rates: it shows the
# rho-epsilon confounding (Stadler 2009; Cusimano & Renner 2010 Syst. Biol.
# 59:458; Hoehna, Stadler, Ronquist & Britton 2011 Mol. Biol. Evol. 28:2577).
# When lambda rho / rho' < r the identity cannot hold, mu is forced to 0, and
# the likelihood drops.  The fits are nested: in the parameters of the fit at
# rho, the fit at rho' > rho is the same model restricted to
# epsilon >= 1 - rho / rho' (its mu' >= 0), a half-space whose edge lies
# inside (0, 1).  The test of that restriction (H0: the pruned tree was
# sampled with rho') has, at the least favourable point epsilon =
# 1 - rho / rho', the 50:50 mixture of chi2(0) and chi2(1) as its null
# distribution (Chernoff 1954; Self & Liang 1987), 5% point qchisq(0.90, 1).
# The rebuild classifies each pair from these outputs and verify_atlas.py
# checks the nesting (boundary iff epsilon at rho < 1 - rho / rho').
#
# The family fit below is reported beside the uniform fit, not compared with
# it: the two likelihoods carry different probabilities of sampling the same
# tips (prod f_i against rho^n), terms that do not depend on the rates, and
# the family fractions are computed from the tips themselves; their
# difference is not a model comparison.
#
# NON-UNIFORM SAMPLING.  rho = tips/described assumes every species is in the
# tree with the same probability.  The deposit contradicts that (the
# per-family molecular fractions printed by the rebuild).  The fit
# 'molecular / family' therefore gives each tip the sampling fraction of its
# Catalogue of Life family (accepted names among the molecular tips of that
# family over the family's accepted species in COL26.9,
# Data/standardized/col269_family_species_counts.csv, capped at 1; a tip with
# no resolved family takes rho = tips/described), with the clade-specific
# likelihood of FitzJohn, Maddison & Otto 2009 (Syst. Biol. 58:595; the
# sampling scheme BAMM and ClaDS use): along a branch, with u = 1 - E,
#     u(t) = r u0 e^{rt} / (r + lambda u0 (e^{rt} - 1))
#     log D(t) = log D0 + r t - 2 log((r + lambda u0 (e^{rt} - 1)) / r),
# D0 = u0 = f at a tip, D = lambda D_left D_right at a node, u taken from the
# first child, and the crown conditioned as above.  With every f equal it
# reproduces the uniform likelihood (checked on every tree, column
# family_uniform_check).  The intervals are
# conditional on the single tree fitted: they exclude topological and dating
# uncertainty, and uniform sampling is an assumption the molecular-pruned
# trees violate (the per-family molecular fraction is reported by the
# rebuild).  Implemented here because TreePar is not installed; checked
# against phytools::fit.bd on every fitted tree of up to CROSSCHECK_MAX_TIPS
# tips (reported in the output).
#
# POLYTOMIES (round F1, finding F1-14).  Consensus summary trees carry nodes
# with more than two children; zero-length internal branches are the same thing
# written as a binary tree.  Every fit records, on its tree after
# di2multi(tol = ZERO_BRANCH_TOL), the internal nodes with more than two
# children, the branching events at them (k - 1 per node) and all branching
# events.  Sensitivity (worked_example_bd_polytomy.csv): each uniform-sampling
# fit is repeated with the polytomies resolved, by keeping one event at the
# node's age and drawing its other k - 2 events uniformly between that age and
# the age of the node's oldest child, for seeds POLY_SEEDS.  Under uniform
# sampling the likelihood depends on the branching times only, so no topology
# needs to be drawn.
#
# CHILD ORDER OF THE FAMILY FIT (finding F1-16).  The clade-specific
# likelihood takes u at a node from one child, so its value depends on the
# order in which children are stored.  The family fit is repeated with u taken
# from the last child instead of the first (worked_example_bd_family_order.csv).
#
# OUTPUT  Paper/inputs/worked_example_bd.csv            one row per fit
#         Paper/inputs/worked_example_bd_selection.csv  the eligible partitions
#         Paper/inputs/worked_example_bd_crosscheck.csv phytools agreement
#         Paper/inputs/worked_example_bd_polytomy.csv   polytomies resolved
#         Paper/inputs/worked_example_bd_family_order.csv  family fit, u from
#                                                       the last child
#
# USAGE   Rscript Paper/code/analysis/worked_example_bd.R [REPO_ROOT]

suppressPackageStartupMessages({ library(ape) })

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  normalizePath(file.path(dirname(f), "..", "..", ".."))
}
P <- function(...) file.path(root, ...)

CUTOFF <- qchisq(0.95, 1) / 2      # profile-likelihood 95% interval
EPS_BOUNDARY <- 1e-6               # an epsilon MLE below this is on the boundary
EPS_MAX <- 1 - 1e-9                # epsilon < 1 (mu < lambda)
# phytools::fit.bd (the cross-check) exhausts 16 GB above ~1e4 tips; it is run
# on every fitted tree up to this size, which validates the likelihood code.
CROSSCHECK_MAX_TIPS <- 15000
# A branch shorter than this (Myr, after rescaling) is a zero-length branch:
# its two ends are one polytomy (ape::di2multi tolerance).
ZERO_BRANCH_TOL <- 1e-6
# Seeds of the polytomy-resolution sensitivity (F1-14).
POLY_SEEDS <- c(1L, 2L, 3L)

basis <- read.csv(P("Data", "coverage_basis.csv"), stringsAsFactors = FALSE)
moltips <- read.csv(P("Data", "molecular_tips.csv"), stringsAsFactors = FALSE,
                    colClasses = "character")
ages <- read.csv(P("Paper", "inputs", "crown_ages.csv"), stringsAsFactors = FALSE)
elig <- subset(basis, is_prokaryotic == "no" & is_chronogram == "yes" &
                 time_scale == "absolute" & partition %in% moltips$partition)
elig$n_molecular <- sapply(elig$partition, function(p)
  sum(moltips$partition == p))
elig <- elig[elig$n_molecular < elig$labels_unique, ]
elig <- elig[order(elig$partition), ]
elig$crown_age_ma <- ages$crown_age_ma[match(elig$partition, ages$partition)]
elig$geometry_root_ma <- ages$geometry_root_ma[match(elig$partition,
                                                     ages$partition)]
elig$time_scale_factor <- 1   # D42: the archived branch lengths, unscaled
elig$rho_C0 <- elig$c0_units / elig$described
elig$rho_full <- elig$labels_unique / elig$described

# ---- family sampling fractions (fit 'molecular / family') -------------------
RESOLVED <- c("accepted", "provisional", "synonym", "accepted_orthographic",
              "provisional_orthographic", "synonym_orthographic")
colres <- read.csv(P("Data", "standardized", "col_resolution.csv"),
                   stringsAsFactors = FALSE, colClasses = "character")
colres <- colres[colres$partition %in% elig$partition &
                 colres$match_type %in% RESOLVED & nzchar(colres$col_family), ]
famcount <- read.csv(P("Data", "standardized", "col269_family_species_counts.csv"),
                     stringsAsFactors = FALSE)
famcount$key <- paste(famcount$kingdom, famcount$family, sep = "|")
family_fractions <- function(part, ids) {
  cr <- colres[colres$partition == part & colres$standardized_id %in% ids, ]
  cr$key <- paste(cr$col_kingdom, cr$col_family, sep = "|")
  agg <- aggregate(list(n_molecular_tips = cr$standardized_id),
                   list(key = cr$key), length)
  nm <- aggregate(list(n_accepted_names = cr$col_id), list(key = cr$key),
                  function(v) length(unique(v)))
  agg <- merge(agg, nm, by = "key")
  agg$n_col_species <- famcount$n_accepted_species[match(agg$key, famcount$key)]
  agg$f <- ifelse(is.na(agg$n_col_species), NA,
                  pmin(1, agg$n_accepted_names / agg$n_col_species))
  agg$partition <- part
  agg$kingdom <- sub("\\|.*", "", agg$key); agg$family <- sub(".*\\|", "", agg$key)
  list(table = agg, tip_key = setNames(cr$key, cr$standardized_id))
}

# ---- likelihood -----------------------------------------------------------
branching_events <- function(tr) {
  n <- Ntip(tr)
  depth <- node.depth.edgelength(tr)
  h <- max(depth[1:n])
  ages <- h - depth[(n + 1):(n + tr$Nnode)]
  kids <- tabulate(tr$edge[, 1], nbins = n + tr$Nnode)[(n + 1):(n + tr$Nnode)]
  ev <- rep(ages, kids - 1)
  sort(pmax(ev, 0), decreasing = TRUE)
}

polytomy_stats <- function(tr) {
  t <- di2multi(tr, tol = ZERO_BRANCH_TOL)
  n <- Ntip(t)
  kids <- tabulate(t$edge[, 1], nbins = n + t$Nnode)[(n + 1):(n + t$Nnode)]
  # the crown node with k children carries k - 1 events, as every other node
  c(n_events = sum(kids - 1L), n_polytomy_nodes = sum(kids > 2L),
    n_polytomy_events = sum((kids - 1L)[kids > 2L]),
    max_children = max(kids))
}

# Branching times with every polytomy resolved at random (F1-14): one event at
# the node's age, its other k - 2 drawn uniformly between that age and the
# age of the node's oldest child.
resolved_events <- function(tr, seed) {
  t <- di2multi(tr, tol = ZERO_BRANCH_TOL)
  n <- Ntip(t); N <- n + t$Nnode
  depth <- node.depth.edgelength(t)
  age <- max(depth[1:n]) - depth
  age[age < 0] <- 0; age[1:n] <- 0
  kids <- tabulate(t$edge[, 1], nbins = N)
  oldest <- numeric(N)
  oc <- tapply(age[t$edge[, 2]], t$edge[, 1], max)
  oldest[as.integer(names(oc))] <- oc
  set.seed(seed)
  inner <- (n + 1L):N
  ex <- rep(inner, pmax(kids[inner] - 2L, 0L))
  ev <- c(age[inner], runif(length(ex), oldest[ex], age[ex]))
  sort(pmax(ev, 0), decreasing = TRUE)
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

# clade-specific sampling: tree t, per-tip sampling fractions f (named by tip
# label).  Vectorised over the nodes of one height level at a time.
prep_tree <- function(t, f, child = "first") {
  t <- reorder(t, "postorder")
  n <- Ntip(t); N <- n + t$Nnode
  depth <- node.depth.edgelength(t)
  age <- max(depth[1:n]) - depth
  age[age < 0] <- 0
  age[1:n] <- 0     # tips are the present, as in branching_events()
  lev <- integer(N)
  for (e in seq_len(nrow(t$edge))) {
    pa <- t$edge[e, 1]; ch <- t$edge[e, 2]
    lev[pa] <- max(lev[pa], lev[ch] + 1L)
  }
  first <- if (child == "first") !duplicated(t$edge[, 1]) else
    !duplicated(t$edge[, 1], fromLast = TRUE)  # the child u is taken from
  elev <- lev[t$edge[, 1]]
  levels <- lapply(split(seq_len(nrow(t$edge)), elev), function(sel) {
    pa <- t$edge[sel, 1]
    upa <- unique(pa)
    list(ch = t$edge[sel, 2], pa_idx = match(pa, upa), upa = upa,
         dt = age[pa] - age[t$edge[sel, 2]],
         first_pos = which(first[sel]), first_pa = pa[first[sel]])
  })
  nk <- tabulate(t$edge[, 1], nbins = N)
  internal <- (n + 1L):N
  extra <- nk[internal] - 1L
  extra[internal == n + 1L] <- nk[n + 1L] - 2L
  list(n = n, N = N, levels = levels, f = unname(f[t$tip.label]),
       root = n + 1L, n_extra = sum(extra))
}

loglik_clade <- function(r, eps, pt) {
  lam <- r / (1 - eps)
  logD <- numeric(pt$N); u <- numeric(pt$N)
  logD[1:pt$n] <- log(pt$f); u[1:pt$n] <- pt$f
  for (lv in pt$levels) {
    ch <- lv$ch
    em <- exp(-r * lv$dt)                    # stable form: divide by e^{r dt}
    den <- r * em + lam * u[ch] * (1 - em)
    utop <- r * u[ch] / den
    ltop <- logD[ch] - r * lv$dt - 2 * log(den / r)
    logD[lv$upa] <- rowsum(ltop, lv$pa_idx)[, 1]
    u[lv$first_pa] <- utop[lv$first_pos]
  }
  logD[pt$root] + pt$n_extra * log(lam) - 2 * log(u[pt$root])
}

fit_clade <- function(pt, x) {
  nll <- function(p) -loglik_clade(exp(p[1]), plogis(p[2]) * EPS_MAX, pt)
  best <- NULL
  for (e0 in c(-4, 0, 2)) {
    r0 <- log((length(x) + 1) / x[1]) / x[1]
    o <- optim(c(log(max(r0, 1e-4)), e0), nll, method = "Nelder-Mead",
               control = list(maxit = 4000, reltol = 1e-12))
    if (is.null(best) || o$value < best$value) best <- o
  }
  list(r = exp(best$par[1]), eps = plogis(best$par[2]) * EPS_MAX, ll = -best$value)
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
  r <- exp(best$par[1]); eps <- plogis(best$par[2]) * EPS_MAX
  list(r = r, eps = eps, ll = -best$value)
}

profile_ci <- function(x, rho, m, which, ll = function(r, e) loglik(r, e, x, rho)) {
  # max over the other parameter at a fixed value of `which`
  prof <- function(v) {
    if (which == "r") {
      f <- function(q) -ll(v, plogis(q) * EPS_MAX)
      o <- optimize(f, c(-25, 25))
    } else {
      f <- function(q) -ll(exp(q), v)
      o <- optimize(f, c(log(1e-8), log(10)))
    }
    -o$objective
  }
  at <- if (which == "r") m$r else m$eps
  g <- function(v) m$ll - prof(v) - CUTOFF
  lo <- hi <- NA
  if (which == "r") {
    lo <- tryCatch(uniroot(g, c(at * 1e-3, at))$root, error = function(e) NA)
    hi <- tryCatch(uniroot(g, c(at, at * 1e3))$root, error = function(e) NA)
  } else {
    lo <- if (at < EPS_BOUNDARY || g(0) <= 0) 0 else
      tryCatch(uniroot(g, c(0, at))$root, error = function(e) NA)
    hi <- if (g(EPS_MAX) <= 0) 1 else
      tryCatch(uniroot(g, c(at, EPS_MAX))$root, error = function(e) NA)
  }
  c(lo, hi)
}

# ---- fits -----------------------------------------------------------------
FITS <- list(c("full", "tips"), c("molecular", "tips"),
             c("molecular", "full_misapplied"), c("molecular", "family"))
rows <- list()
check <- NULL
famtabs <- list()
poly <- list(); ford <- list()
for (i in seq_len(nrow(elig))) {
  b <- elig[i, ]
  tr <- read.tree(P("Data", "standardized", "trees", b$canonical_tree))
  tr$edge.length <- tr$edge.length * b$time_scale_factor
  mol <- intersect(moltips$standardized_id[moltips$partition == b$partition],
                   tr$tip.label)
  trees <- list(full = tr, molecular = keep.tip(tr, mol))
  for (fc in FITS) {
    tn <- fc[1]; conv <- fc[2]
    t <- trees[[tn]]
    x <- branching_events(t)
    rho_t <- min(Ntip(t) / b$described, 1)
    uni_check <- NA_real_; n_tip_family <- NA_integer_
    if (conv == "family") {
      ff <- family_fractions(b$partition, t$tip.label)
      famtabs[[b$partition]] <- ff$table
      fk <- ff$table$f[match(ff$tip_key[t$tip.label], ff$table$key)]
      n_tip_family <- sum(!is.na(fk))
      fk[is.na(fk)] <- rho_t
      pt <- prep_tree(t, setNames(fk, t$tip.label))
      # with every f equal to rho the clade likelihood is the uniform one
      pu <- pt; pu$f[] <- rho_t
      m0 <- fit(x, rho_t)
      uni_check <- abs(loglik_clade(m0$r, m0$eps, pu) - m0$ll)
      m <- fit_clade(pt, x)
      llf <- function(r, e) loglik_clade(r, e, pt)
      ci_r <- profile_ci(x, NA, m, "r", llf)
      ci_e <- profile_ci(x, NA, m, "eps", llf)
      rho <- NA_real_
    } else {
      rho <- if (conv == "tips") rho_t else min(Ntip(trees$full) / b$described, 1)
      m <- fit(x, rho)
      ci_r <- profile_ci(x, rho, m, "r")
      ci_e <- profile_ci(x, rho, m, "eps")
    }
    lam <- m$r / (1 - m$eps)
    ps <- polytomy_stats(t)
    if (conv == "family") {
      po <- prep_tree(t, setNames(fk, t$tip.label), child = "last")
      mo <- fit_clade(po, x)
      ford[[length(ford) + 1]] <- data.frame(
        partition = b$partition, tree = tn, rho_convention = conv,
        child = c("first", "last"),
        r = signif(c(m$r, mo$r), 6), epsilon = signif(c(m$eps, mo$eps), 6),
        loglik = round(c(m$ll, mo$ll), 4),
        loglik_last_at_first_mle = round(loglik_clade(m$r, m$eps, po), 4),
        stringsAsFactors = FALSE)
    } else {
      for (sd in POLY_SEEDS) {
        xs <- resolved_events(t, sd)
        ms <- fit(xs, rho)
        lams <- ms$r / (1 - ms$eps)
        poly[[length(poly) + 1]] <- data.frame(
          partition = b$partition, tree = tn, rho_convention = conv,
          seed = sd, rho = signif(rho, 8), lambda = signif(lams, 6),
          mu = signif(ms$eps * lams, 6), r = signif(ms$r, 6),
          epsilon = signif(ms$eps, 6), loglik = round(ms$ll, 4),
          stringsAsFactors = FALSE)
      }
    }
    rows[[length(rows) + 1]] <- data.frame(
      partition = b$partition, tree = tn,
      canonical_tree = b$canonical_tree, n_tips = Ntip(t),
      crown_age_ma = round(x[1], 4),
      time_scale_factor = signif(b$time_scale_factor, 8),
      rho_convention = conv, rho = signif(rho, 8),
      lambda = signif(lam, 6), mu = signif(m$eps * lam, 6),
      r = signif(m$r, 6), r_lo = signif(ci_r[1], 6), r_hi = signif(ci_r[2], 6),
      epsilon = signif(m$eps, 6), epsilon_lo = signif(ci_e[1], 6),
      epsilon_hi = signif(ci_e[2], 6),
      epsilon_ci = if (m$eps < EPS_BOUNDARY) "estimate at 0" else "interior",
      loglik = round(m$ll, 4),
      n_tips_family_fraction = n_tip_family,
      family_uniform_check = signif(uni_check, 3),
      n_events = ps[["n_events"]], n_polytomy_nodes = ps[["n_polytomy_nodes"]],
      n_polytomy_events = ps[["n_polytomy_events"]],
      max_children = ps[["max_children"]],
      stringsAsFactors = FALSE)
    cat(sprintf("%-12s %-9s %-13s n=%7d rho=%.5f r=%.5f [%.5f, %.5f] eps=%.4f [%.4f, %.4f]\n",
                b$partition, tn, conv, Ntip(t), rho, m$r, ci_r[1], ci_r[2],
                m$eps, ci_e[1], ci_e[2]))
    if (conv != "family" && Ntip(t) <= CROSSCHECK_MAX_TIPS &&
        requireNamespace("phytools", quietly = TRUE)) {
      # cross-check of the likelihood against phytools::fit.bd (same model,
      # same conditioning): maximum-likelihood rates must agree.
      # phytools tests ultrametricity at a tolerance tighter than the
      # rounding of the deposited branch lengths (root-to-tip CV ~1e-7);
      # 'extend' lengthens terminal branches by that rounding only
      tu <- phytools::force.ultrametric(multi2di(t, random = FALSE),
                                        method = "extend")
      pf <- phytools::fit.bd(tu, rho = rho)
      check <- rbind(check, data.frame(
        partition = b$partition, tree = tn, rho_convention = conv,
        rho = signif(rho, 8),
        this_lambda = signif(lam, 6), this_mu = signif(m$eps * lam, 6),
        phytools_lambda = signif(pf$b, 6), phytools_mu = signif(pf$d, 6),
        stringsAsFactors = FALSE))
      print(tail(check, 1))
    }
  }
}
out <- do.call(rbind, rows)
dir.create(P("Paper", "inputs"), showWarnings = FALSE)
write.csv(out, P("Paper", "inputs", "worked_example_bd.csv"), row.names = FALSE, na = "")
sel <- elig[, c("partition", "canonical_tree", "described", "c0_units",
                "labels_unique", "n_molecular", "crown_age_ma",
                "geometry_root_ma", "time_scale_factor", "rho_C0", "rho_full")]
write.csv(sel, P("Paper", "inputs", "worked_example_bd_selection.csv"),
          row.names = FALSE)
fam <- do.call(rbind, famtabs)
fam <- fam[order(fam$partition, fam$kingdom, fam$family),
           c("partition", "kingdom", "family", "n_molecular_tips",
             "n_accepted_names", "n_col_species", "f")]
fam$f <- signif(fam$f, 8)
write.csv(fam, P("Paper", "inputs", "worked_example_bd_family_fractions.csv"),
          row.names = FALSE, na = "")
write.csv(do.call(rbind, poly), P("Paper", "inputs", "worked_example_bd_polytomy.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, ford),
          P("Paper", "inputs", "worked_example_bd_family_order.csv"),
          row.names = FALSE)
if (!is.null(check))
  write.csv(check, P("Paper", "inputs", "worked_example_bd_crosscheck.csv"),
            row.names = FALSE)
cat("wrote Paper/inputs/worked_example_bd*.csv\n")
