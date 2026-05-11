.atlas_base <- function() {
  getOption(
    "phyloatlas.base_url",
    "https://raw.githubusercontent.com/franciscorichter/phylo-species-atlas/main"
  )
}

.cache <- new.env(parent = emptyenv())

.load_dictionary <- function() {
  if (is.null(.cache$dictionary)) {
    .cache$dictionary <- utils::read.csv(
      paste0(.atlas_base(), "/standardized/dictionary.csv"),
      stringsAsFactors = FALSE
    )
  }
  .cache$dictionary
}

.load_metadata <- function() {
  if (is.null(.cache$metadata)) {
    md <- utils::read.csv(
      paste0(.atlas_base(), "/standardized/metadata.csv"),
      stringsAsFactors = FALSE
    )
    prov <- utils::read.csv(
      paste0(.atlas_base(), "/data_provenance.csv"),
      stringsAsFactors = FALSE
    )
    md$name <- sub("\\.nwk$", "", md$filename)
    .cache$metadata <- md
    .cache$provenance <- prov
  }
  .cache$metadata
}

.load_provenance <- function() {
  .load_metadata()
  .cache$provenance
}

#' Clear the in-memory cache
#'
#' The dictionary and metadata files are downloaded once per R session and
#' cached. Call this if you want to force a re-download (for example after
#' switching `phyloatlas.base_url`).
#'
#' @return Invisibly `NULL`.
#' @export
atlas_clear_cache <- function() {
  rm(list = ls(.cache), envir = .cache)
  invisible(NULL)
}
