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
phyfromnewick(file = "", text)
phyfromphylo(phylo)
}
\arguments{
  \item{file}{
a textfile containing the newick phylogeny, which must end with a semi-colon (\code{phyfromnewick} only)}
  \item{text}{
a character variable containing the newick phylogeny, which must end with a semi-colon. Exactly one of these arguments should be supplied (\code{phyfromnewick} only)}
  \item{phylo}{an R-variable containing a "phylo" object, as used, for example, in the \code{ape} package (\code{phyfromphylo} only).}
}
\details{
 Before passing the text to \code{read.newick}, I strip out all the spaces, to avoid problems.
 	
\code{phyfromnewick} contains a function \code{read.newick} which "reads Newick style tree with branch lengths into memory as an \code{ape} "phylo" object", and was written by Liam J. Revell in 2011, and downloaded on 24th January 2014. See references for url.

 Examples are given under \code{\link{phyreg}}.
}
\value{
\item{phy}{The phylogenetic vector in internal format, with an element for every node, except the root, but including species. If \code{jj==phy[[ii]}, then \code{jj} is the parent-node of \code{ii}. It is a requirement for the internal format that \code{jj>ii} in every case (it is \emph{not} true of "phylo" objects that a parent's ID has to be greater than its offspring's). Suitable for supplying as the \code{phydata} argument of \code{phyreg}}
\item{hts}{The heights of each node, a vector with an element for every node. Thus it is longer by one than $phy. Suitable for supplying as the \code{heightsdata} argument of \code{phyreg}}
\item{orphylo}{The "phylo" style object created by \code{read.newick}. (\code{phylofromnewick} only)}
}
\references{
\url{https://github.com/liamrevell/phytools/blob/master/R/read.newick.R} %% I got it from somewhere else
%% https://github.com/cboettig/Comparative-Phylogenetics/blob/master/RevellExtensions/R/read.newick.R, but it disappeared by 2018_04_05. I haven't used the updated version, I'm afraid
}
\author{
Alan Grafen,
with an internal function (\code{read.newick}) written by Liam J. Revell in 2011, downloaded on 24th January 2014. See references for url.
}
\seealso{
\code{\link{phyreg}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Conversion }
\keyword{ Phylogeny }% __ONLY ONE__ keyword per line
