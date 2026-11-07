.atlas_base <- function() {
  getOption(
    "phyloatlas.base_url",
    "https://raw.githubusercontent.com/franciscorichter/phylo-species-atlas/main"
  )
}

.cache <- new.env(parent = emptyenv())

# Read a CSV/Newick resource from a URL with an enforced session-local
# timeout. Returns the read value on success, or NULL with a single
# diagnostic message() on failure. Used to satisfy the CRAN policy that
# packages must "fail gracefully if a resource is not available" without
# throwing or hanging CRAN's check infrastructure.
.try_url <- function(url, reader, what) {
  old <- options(timeout = max(30, getOption("timeout", 60L)))
  on.exit(options(old), add = TRUE)
  tryCatch(
    reader(url),
    error = function(e) {
      message(sprintf(
        "phyloatlas: could not fetch %s from %s. Underlying error: %s",
        what, url, conditionMessage(e)
      ))
      NULL
    }
  )
}

.load_dictionary <- function() {
  if (is.null(.cache$dictionary)) {
    dict <- .try_url(
      paste0(.atlas_base(), "/standardized/dictionary.csv"),
      function(u) utils::read.csv(u, stringsAsFactors = FALSE),
      "dictionary"
    )
    if (is.null(dict)) return(NULL)
    if (anyDuplicated(dict$id)) {
      warning(
        "Atlas dictionary contains duplicate ids; using first match for each.",
        call. = FALSE
      )
    }
    .cache$dictionary <- dict
  }
  .cache$dictionary
}

.load_metadata <- function() {
  if (is.null(.cache$metadata)) {
    md <- .try_url(
      paste0(.atlas_base(), "/standardized/metadata.csv"),
      function(u) utils::read.csv(u, stringsAsFactors = FALSE),
      "tree metadata"
    )
    prov <- .try_url(
      paste0(.atlas_base(), "/data_provenance.csv"),
      function(u) utils::read.csv(u, stringsAsFactors = FALSE),
      "data provenance"
    )
    if (is.null(md) || is.null(prov)) return(NULL)
    md$name <- sub("\\.nwk$", "", md$filename)
    .cache$metadata <- md
    .cache$provenance <- prov
  }
  .cache$metadata
}

.load_provenance <- function() {
  if (is.null(.load_metadata())) return(NULL)
  .cache$provenance
}

# The trees the atlas lists but does not redistribute, because their sources
# state no licence: registers/withheld_trees.csv under the base URL (data
# release v1.1.3 and later), else the copy bundled with the package. Returns a
# data frame with columns name, study, source_url, licence (possibly empty).
.load_withheld <- function() {
  if (is.null(.cache$withheld)) {
    url <- paste0(.atlas_base(), "/registers/withheld_trees.csv")
    old <- options(timeout = max(30, getOption("timeout", 60L)))
    on.exit(options(old), add = TRUE)
    reg <- tryCatch(
      suppressWarnings(utils::read.csv(url, stringsAsFactors = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(reg) && all(c("release_path", "source_url") %in% names(reg))) {
      trees <- startsWith(reg$release_path, "standardized/trees/")
      reg <- reg[trees, , drop = FALSE]
      out <- data.frame(
        name = sub("\\.nwk$", "", basename(reg$release_path)),
        study = if ("study" %in% names(reg)) reg$study else NA_character_,
        source_url = reg$source_url,
        licence = if ("licence" %in% names(reg)) reg$licence else NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      path <- system.file("extdata", "withheld_trees.csv", package = "phyloatlas")
      out <- if (nzchar(path)) {
        utils::read.csv(path, stringsAsFactors = FALSE)
      } else {
        data.frame(name = character(0), study = character(0),
                   source_url = character(0), licence = character(0))
      }
    }
    .cache$withheld <- out
  }
  .cache$withheld
}

#' Clear the in-memory cache
#'
#' The dictionary and metadata files are downloaded once per R session and
#' cached. Call this if you want to force a re-download (for example after
#' switching `phyloatlas.base_url`).
#'
#' @return Invisibly `NULL`. Called for its side effect of emptying the
#'   internal cache environment.
#'
#' @family atlas
#' @examples
#' # Safe to run unconditionally — just empties the in-session cache.
#' atlas_clear_cache()
#' @export
atlas_clear_cache <- function() {
  rm(list = ls(.cache), envir = .cache)
  invisible(NULL)
}
