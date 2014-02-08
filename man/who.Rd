\name{Identifying higher nodes}
\alias{who}

\title{
Identifying higher nodes in an internal-format phylogeny}
\description{
The internal ID numbers of species are just their position in the main datafile, which should be the same as the order in internal, taxonomic vector and newick formats for phylogenies. Identifying higher nodes from their ID numbers can be more difficult. \code{who} allows you to name a node, or list of nodes, and a phylogeny, and it provides information to help you identify the node(s). This might be useful when looking at short regression output, or the plot from \code{phyreg}.
}
\usage{
who(phy, nodes, find = "mm")
}
\arguments{
  \item{phy}{
This is an internal-format phylogeny, as supplied in \code{phyreg} output as \code{$fullphy} for the full phylogeny or \code{$usedphy} for the reduced phylogeny (reflecting species omitted either because of the \code{subset} argument of \code{phyreg} or because of missing values). The full phylogeny will usually be the right one to use, because all the usual output uses those originalIDs, rather than the reduced ones.
}
  \item{nodes}{
an integer representing an node ID number, or a list, for example \code{c(101,123)} or \code{101:137}. If not specified, it will be taken to be the list of all higher nodes.
}
  \item{find}{if \code{"sp"} then all the species below the node are printed, if \code{"mm"} then the minimum and maximum of those species numbers are printed, and if \code{"dg"} (or anything else) then the IDs of the daughter nodes are printed.
}
}
\value{
	If \code{find="sp"} or \code{find="dg"}, the results are printed out, and nothing is returned.
	
	If \code{find="mm"}, a data-frame is returned with three columns: \code{ID} for the number of the higher node, \code{minsp} for the first species, and \code{maxsp} for the last species.
	
}
\details{
	If the species below a higher node are a contiguous set, then \code{find="mm"} (for min and max) will conveniently tell you the first and last species in that set. If your higher nodes have interdigitating sets of species, then you may need to choose \code{"sp"} to see the complete set. For drawing out your own phylogeny, \code{find="dg"} for the daughter nodes could be helpful.
	
	Examples of use are provided in \code{\link{phyreg}}
	}
\author{
Alan Grafen
}

\seealso{
\code{\link{phyreg}} uses these phylogenies and produces the output that \code{who} can help to interpret. 
}
