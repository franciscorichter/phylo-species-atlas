#!/usr/bin/env Rscript
# Extract crown-age 95% intervals from source trees for Figure 2 error bars.
#
# Strategy (avoids fragile treeio dependencies):
#   beast_mcc_annotation -> regex over the file: find all "height_95%_HPD={low,high},height=H"
#                           the max-H node is the root; its HPD is the crown-age HPD.
#   posterior_sample / bootstrap_treepl -> read each "tree NAME = ... ;" line, strip BEAST
#                           [&...] annotations and rooting markers, parse the de-annotated
#                           Newick with ape::read.tree, compute crown age per tree,
#                           report 2.5%/50%/97.5% quantiles across replicates.
#
# Output: Paper/inputs/fig2_crown_age_hpds.csv

suppressMessages(library(ape))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "."
script_dir <- normalizePath(dirname(script_path), mustWork = FALSE)
ATLAS_DATA   <- normalizePath(file.path(script_dir, "..", "..", "..", "Data"), mustWork = FALSE)
PAPER_INPUTS <- normalizePath(file.path(script_dir, "..", "..", "inputs"), mustWork = FALSE)
INV_FILE   <- file.path(PAPER_INPUTS, "fig2_hpd_inventory.csv")
OUT_FILE   <- file.path(PAPER_INPUTS, "fig2_crown_age_hpds.csv")

inv <- read.csv(INV_FILE, stringsAsFactors = FALSE)
inv <- inv[inv$hpd_method != "none" & nzchar(inv$source_file), ]

# ---- BEAST MCC: regex-extract root HPD --------------------------------------
extract_beast_mcc <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # Match: height_95%_HPD={low,high},height=H   (height may precede HPD; allow either order)
  # Two patterns:
  m1 <- gregexpr(
    "height_95%_HPD=\\{([0-9eE.+-]+),([0-9eE.+-]+)\\}[^]]*?height=([0-9eE.+-]+)",
    txt, perl = TRUE)
  m2 <- gregexpr(
    "height=([0-9eE.+-]+)[^]]*?height_95%_HPD=\\{([0-9eE.+-]+),([0-9eE.+-]+)\\}",
    txt, perl = TRUE)

  parse_caps <- function(matches, hpd_first) {
    if (matches[[1]][1] == -1) return(NULL)
    starts <- matches[[1]]
    lengths <- attr(matches[[1]], "match.length")
    out <- lapply(seq_along(starts), function(k) {
      sub <- substr(txt, starts[k], starts[k] + lengths[k] - 1)
      lo <- as.numeric(sub("^.*?height_95%_HPD=\\{([0-9eE.+-]+),.*$", "\\1", sub, perl = TRUE))
      hi <- as.numeric(sub("^.*?height_95%_HPD=\\{[0-9eE.+-]+,([0-9eE.+-]+)\\}.*$", "\\1", sub, perl = TRUE))
      h  <- as.numeric(sub("^.*?height=([0-9eE.+-]+).*$", "\\1", sub, perl = TRUE))
      c(h, lo, hi)
    })
    do.call(rbind, out)
  }

  a <- parse_caps(m1, TRUE)
  b <- parse_caps(m2, FALSE)
  cand <- rbind(a, b)
  if (is.null(cand) || nrow(cand) == 0) {
    stop("no height_95%_HPD annotations found")
  }
  cand <- unique(cand)
  root <- cand[which.max(cand[, 1]), , drop = TRUE]
  list(point = root[1], low = root[2], high = root[3], n = 1L,
       n_annotated_nodes = nrow(cand), err = NA_character_)
}

# ---- multi-tree file: strip & parse each Newick ----------------------------
strip_beast_anno <- function(s) {
  # remove [& ... ] blocks (possibly nested-ish but BEAST doesn't nest them within a single annotation)
  s <- gsub("\\[&[^]]*\\]", "", s, perl = TRUE)
  s
}

extract_multitree <- function(path) {
  # Strategy: read the whole file, strip every [&...] annotation block,
  # write the cleaned text to a temp file, then let ape::read.nexus or
  # ape::read.tree do the heavy lifting (translate blocks, multiple trees).
  raw <- readLines(path, warn = FALSE)
  cleaned <- gsub("\\[&[^]]*\\]", "", raw, perl = TRUE)
  tmp <- tempfile(fileext = ".tre")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(cleaned, tmp)

  is_nexus <- any(grepl("#NEXUS", cleaned[1:min(5, length(cleaned))], ignore.case = TRUE))
  trs <- if (is_nexus) {
    tryCatch(read.nexus(tmp), error = function(e) NULL)
  } else {
    tryCatch(read.tree(tmp), error = function(e) NULL)
  }
  # Fallbacks
  if (is.null(trs)) {
    trs <- tryCatch(read.tree(tmp), error = function(e) NULL)
  }
  if (is.null(trs)) {
    trs <- tryCatch(read.nexus(tmp), error = function(e) NULL)
  }
  if (is.null(trs)) stop("failed to parse cleaned tree file")
  if (inherits(trs, "phylo")) trs <- list(trs)
  ca <- sapply(trs, function(tr) {
    tryCatch(max(branching.times(tr)), error = function(e) NA_real_)
  })
  ca <- ca[!is.na(ca) & ca > 0]
  if (length(ca) == 0) stop("no parseable crown ages")
  q <- quantile(ca, probs = c(0.025, 0.5, 0.975), names = FALSE)
  list(point = q[2], low = q[1], high = q[3], n = length(ca),
       n_annotated_nodes = NA_integer_, err = NA_character_)
}

out <- data.frame()
for (i in seq_len(nrow(inv))) {
  row <- inv[i, ]
  full <- file.path(ATLAS_DATA, row$source_file)
  cat(sprintf("[%d/%d] %s (%s) ... ", i, nrow(inv), row$dataset, row$hpd_method))
  if (!file.exists(full)) { cat("MISSING\n"); next }
  res <- tryCatch({
    if (row$hpd_method == "beast_mcc_annotation") extract_beast_mcc(full)
    else extract_multitree(full)
  }, error = function(e) list(point = NA, low = NA, high = NA, n = 0L,
                              n_annotated_nodes = NA_integer_, err = conditionMessage(e)))
  if (!is.na(res$point)) {
    cat(sprintf("crown=%.2f [%.2f, %.2f] (n=%d)\n", res$point, res$low, res$high, res$n))
  } else {
    cat(sprintf("ERROR: %s\n", res$err))
  }
  out <- rbind(out, data.frame(
    dataset = row$dataset, hpd_method = row$hpd_method,
    crown_age_point = res$point, crown_age_low = res$low, crown_age_high = res$high,
    n_trees = res$n, source_file = row$source_file,
    notes = if (is.na(res$err)) "" else res$err,
    stringsAsFactors = FALSE
  ))
}

write.csv(out, OUT_FILE, row.names = FALSE)
cat("\nWrote ", OUT_FILE, "\n", sep = "")
print(out[, c("dataset", "hpd_method", "crown_age_point", "crown_age_low", "crown_age_high", "n_trees", "notes")])
