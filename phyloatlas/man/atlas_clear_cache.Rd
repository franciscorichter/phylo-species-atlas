% Generated manually; matches R/utils.R
\name{atlas_clear_cache}
\alias{atlas_clear_cache}
\title{Clear the in-memory cache}
\usage{
atlas_clear_cache()
}
\value{
Invisibly \code{NULL}.
}
\description{
The dictionary and metadata files are downloaded once per R session and
cached. Call this if you want to force a re-download (for example after
switching \code{phyloatlas.base_url}).
}
