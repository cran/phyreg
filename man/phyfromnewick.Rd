\name{Translating phylogenies}
\alias{phyfromnewick}
\alias{phyfromphylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using phylogenies formatted in standard but "other" ways
}
\description{
The two functions convert phylogenies in a standard format (newick or phylo) into a form that can be supplied as the argument \code{phydata} of \code{phyreg}. If node heights are specified, these are also provided in the value returned.
}
\usage{
phyfromnewick(file = "", str,datatype="Branch_lengths")
phyfromphylo(phylo)
}
\arguments{
  \item{file}{
a textfile containing the newick phylogeny
}
  \item{str}{
a character variable containing the newick phylogeny. Exactly one of these arguments should be supplied. (\code{phyfromnewick} only)
}
\item{datatype}{The newick object can give branch lengths (by including a colon and a number after each species and after each higher node except the root), or node heights (by doing so for every node including the root). The function looks for the kind of data suggested by this argument. If there are no numbers, or not all required nodes have them, then they are simply discarded and no heights are returned. (\code{phyfromnewick} only)}
  \item{phylo}{an R-name containing a "phylo" object, as used, for example, in the \code{ape} package The phylo object must have \code{edge.length} values if heights are to be provided by \code{phyfromphylo}. (\code{phyfromphylo} only).}
}
\details{
newick and phylo formats are both widely-used standards. Examples of the functions are given under \code{\link{phyreg}}.

Both standards can be constructed to contain singleton higher nodes i.e. a higher node with only one daughter. \code{phyfromphylo} leaves the singletons in. \code{phyfromnewick} strips them out. \code{phyreg} will use either.
}
\value{
\item{phy}{The phylogenetic vector in internal format, with an element for every node, except the root, but including species. If \code{jj==phy[[ii]]}, then \code{jj} is the parent-node of \code{ii}. It is a requirement for the internal format that \code{jj>ii} in every case (it is \emph{not} true of "phylo" objects that a parent's ID has to be greater than its offspring's). Suitable for supplying as the \code{phydata} argument of \code{phyreg}}
\item{hts}{
The heights of each node, a vector with an element for every node. Thus it is longer by one than \code{$phy}. Suitable for supplying as the \code{heightsdata} argument of \code{phyreg}. NOTE: this will be returned only if the information is available in the inputs.
}
\item{firstphy}{A \code{newick} object may have higher nodes with only one daughter. \code{firstphy} contains the corresponding internally formatted phylogeny. \code{phy} itself has all singletons removed, and the numbers of the higher nodes are therefore altered. (\code{phylofromnewick} only)}
\item{sd}{Species data as read from the \code{newick} object. It contains a name (Species_1 etc if none are supplied), the ID number, and the additional data that, if numeric and provided for all relevant nodes, will be interpreted as branch lengths or heights. This may be useful if numerical data is not getting through for some reason (\code{phylofromnewick} only)}
\item{oh}{Higher node data as read from the \code{newick} object. This contains the higher node names (Higher_Node_1 etc by default), their ID numbers \emph{in \code{firstphy} not in \code{phy}}, and the additional data. Again, perhaps useful for checking. (\code{phylofromnewick} only)}
}
\references{ For newick format, see
 \url{http://en.wikipedia.org/wiki/Newick_format}. For phylo objects, so far as I can see invented by Diaconis, P.W. and S.P. Holmes (1998, PNAS 95, 14600--14602), see the \code{ape} package (\url{http://cran.r-project.org/web/packages/ape/}).
}
\author{
Alan Grafen
}
\seealso{
\code{\link{phyreg}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Conversion }
\keyword{ Phylogeny }% __ONLY ONE__ keyword per line
