#' List trees available in the Phylo-Species Atlas
#'
#' Returns a data frame describing every standardized tree in the atlas:
#' name, study, number of tips, whether the tree is time-calibrated, plus
#' provenance fields (year, journal, DOI, coverage) when available.
#'
#' @return A data frame with one row per tree, ordered by group.
#'
#' @examples
#' \dontrun{
#' trees <- list_trees()
#' head(trees)
#' subset(trees, dated & ntips > 1000)
#' }
#' @export
list_trees <- function() {
  md <- .load_metadata()
  prov <- .load_provenance()
  prov_key <- .provenance_key(md$group, prov$group)
  prov_cols <- setdiff(
    names(prov),
    c("group", "tree_name", "tips", "dated", "format_original",
      "species_count_source")
  )
  joined <- prov[match(prov_key, prov$group), prov_cols, drop = FALSE]
  if ("study" %in% names(joined) && "study" %in% names(md)) {
    names(joined)[names(joined) == "study"] <- "study_full"
  }
  out <- cbind(md, joined)
  cols <- c(
    "name", "group", "study", "ntips", "dated",
    "year", "journal", "doi", "crown_ma",
    "described_species", "coverage_pct",
    "data_source", "download_url", "methods_brief", "notes",
    "study_full"
  )
  cols <- intersect(cols, names(out))
  out <- out[, cols, drop = FALSE]
  out <- out[order(out$group, out$name), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Map a metadata group to the canonical provenance group.
# - direct match wins
# - "birds" -> "birds_mctavish" (the standardized tree is McTavish 2025)
# - "condamine_*" -> "condamine"
.provenance_key <- function(md_groups, prov_groups) {
  key <- md_groups
  key[key == "birds"] <- "birds_mctavish"
  is_cond <- startsWith(key, "condamine_")
  key[is_cond] <- "condamine"
  miss <- !(key %in% prov_groups)
  key[miss] <- NA_character_
  key
}
