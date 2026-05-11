% Generated manually; matches R/list_trees.R
\name{list_trees}
\alias{list_trees}
\title{List trees available in the Phylo-Species Atlas}
\usage{
list_trees()
}
\value{
A data frame with one row per tree, ordered by group.
}
\description{
Returns a data frame describing every standardized tree in the atlas:
name, study, number of tips, whether the tree is time-calibrated, plus
provenance fields (year, journal, DOI, coverage) when available.
}
\examples{
\dontrun{
trees <- list_trees()
head(trees)
subset(trees, dated & ntips > 1000)
}
}
