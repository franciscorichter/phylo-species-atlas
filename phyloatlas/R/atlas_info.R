#' Get metadata for a single tree
#'
#' @param name Tree name without the `.nwk` extension (e.g. `"mammals"`,
#'   `"birds"`, `"seed_plants"`). Must be a length-1 non-`NA` character
#'   string. Use [list_trees()] to see all available names.
#'
#' @return A one-row data frame with all metadata columns from
#'   [list_trees()], or `NULL` (with a warning) if no tree by that name
#'   exists. Also returns `NULL` (with a diagnostic message, no warning)
#'   if the atlas metadata cannot be downloaded.
#'
#' @family atlas
#' @examples
#' \donttest{
#' info <- atlas_info("mammals")
#' if (!is.null(info)) info
#' }
#' @export
atlas_info <- function(name) {
  stopifnot(
    is.character(name), length(name) == 1L, !is.na(name), nzchar(name)
  )
  name <- sub("\\.nwk$", "", name)
  trees <- list_trees()
  if (is.null(trees)) return(NULL)
  row <- trees[trees$name == name, , drop = FALSE]
  if (nrow(row) == 0L) {
    warning("No tree named '", name, "'. See list_trees().", call. = FALSE)
    return(NULL)
  }
  rownames(row) <- NULL
  row
}
