#' Load a tree from the Phylo-Species Atlas
#'
#' Downloads a standardized Newick tree from the atlas and resolves its
#' integer tip IDs to species names using the shared dictionary.
#'
#' @param name Tree name without the `.nwk` extension (e.g. `"mammals"`,
#'   `"birds"`, `"seed_plants"`, `"condamine_Vangidae"`). Use [list_trees()]
#'   to see all available names.
#' @param resolve_labels If `TRUE` (default) tip labels are replaced with
#'   standardized species names from `dictionary.csv`. Set to `FALSE` to
#'   keep the raw integer IDs (faster, avoids downloading the dictionary).
#'
#' @return An object of class `"phylo"` from the \pkg{ape} package. If
#'   `resolve_labels = TRUE` (the default), tip labels are species names
#'   from the standardized dictionary; otherwise tip labels are integer
#'   IDs as character strings.

#' @family atlas
#'
#' @examples
#' # Offline demo using a small bundled tree (does not hit the network):
#' demo_path <- system.file("extdata", "tree_demo.nwk", package = "phyloatlas")
#' tree <- ape::read.tree(demo_path)
#' tree
#'
#' \donttest{
#' # Live atlas fetch (requires internet):
#' tree <- try(load_atlas_tree("mammals"), silent = TRUE)
#' if (!inherits(tree, "try-error")) plot(tree, show.tip.label = FALSE)
#'
#' # Keep integer IDs to skip the 18 MB dictionary download:
#' tree <- try(load_atlas_tree("birds", resolve_labels = FALSE), silent = TRUE)
#' }
#' @export
load_atlas_tree <- function(name, resolve_labels = TRUE) {
  stopifnot(
    is.character(name), length(name) == 1L, !is.na(name), nzchar(name)
  )
  name <- sub("\\.nwk$", "", name)
  url <- paste0(.atlas_base(), "/standardized/trees/", name, ".nwk")
  old <- options(timeout = max(30, getOption("timeout", 60L)))
  on.exit(options(old), add = TRUE)
  tree <- tryCatch(
    ape::read.tree(url),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("cannot open URL|HTTP|Timeout|timeout|resolve host", msg)) {
        stop(
          "Could not fetch tree '", name, "' (network resource unavailable). ",
          "Underlying error: ", msg,
          call. = FALSE
        )
      }
      stop(
        "Could not load tree '", name, "'. ",
        "Check the name with list_trees(). Underlying error: ", msg,
        call. = FALSE
      )
    }
  )
  if (resolve_labels) {
    dict <- .load_dictionary()
    if (!is.null(dict)) {
      ids <- suppressWarnings(as.integer(tree$tip.label))
      mapped <- dict$standardized_name[match(ids, dict$id)]
      keep <- !is.na(mapped)
      tree$tip.label[keep] <- mapped[keep]
    }
  }
  tree
}
