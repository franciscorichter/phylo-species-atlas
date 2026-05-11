% Generated manually; matches R/atlas_info.R
\name{atlas_info}
\alias{atlas_info}
\title{Get metadata for a single tree}
\usage{
atlas_info(name)
}
\arguments{
\item{name}{Tree name without the \code{.nwk} extension.}
}
\value{
A one-row data frame with all available metadata fields, or
\code{NULL} (with a warning) if no tree by that name exists.
}
\description{
Get metadata for a single tree
}
\examples{
\dontrun{
atlas_info("mammals")
}
}
