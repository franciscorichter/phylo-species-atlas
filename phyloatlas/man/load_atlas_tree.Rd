% Generated manually; matches R/load_atlas_tree.R
\name{load_atlas_tree}
\alias{load_atlas_tree}
\title{Load a tree from the Phylo-Species Atlas}
\usage{
load_atlas_tree(name, resolve_labels = TRUE)
}
\arguments{
\item{name}{Tree name without the \code{.nwk} extension (e.g.
\code{"mammals"}, \code{"birds"}, \code{"seed_plants"},
\code{"condamine_Vangidae"}). Use \code{\link{list_trees}} to see all
available names.}

\item{resolve_labels}{If \code{TRUE} (default) tip labels are replaced
with standardized species names from \code{dictionary.csv}. Set to
\code{FALSE} to keep the raw integer IDs (faster, avoids downloading
the dictionary).}
}
\value{
An object of class \code{"phylo"} from the \pkg{ape} package.
}
\description{
Downloads a standardized Newick tree from the atlas and resolves its
integer tip IDs to species names using the shared dictionary.
}
\examples{
\dontrun{
tree <- load_atlas_tree("mammals")
plot(tree, show.tip.label = FALSE)

# Keep integer IDs to avoid loading the 18 MB dictionary
tree <- load_atlas_tree("birds", resolve_labels = FALSE)
}
}
