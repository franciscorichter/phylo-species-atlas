#' Get metadata for a single tree
#'
#' @param name Tree name without the `.nwk` extension.
#'
#' @return A one-row data frame with all available metadata fields, or
#'   `NULL` (with a warning) if no tree by that name exists.
#'
#' @examples
#' \donttest{
#' info <- try(atlas_info("mammals"), silent = TRUE)
#' if (!inherits(info, "try-error")) info
#' }
#' @export
atlas_info <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
  name <- sub("\\.nwk$", "", name)
  trees <- list_trees()
  row <- trees[trees$name == name, , drop = FALSE]
  if (nrow(row) == 0L) {
    warning("No tree named '", name, "'. See list_trees().", call. = FALSE)
    return(NULL)
  }
  rownames(row) <- NULL
  row
}
