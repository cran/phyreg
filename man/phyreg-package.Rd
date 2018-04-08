\name{phyreg-package}
\alias{phyreg-package}
\docType{package}
\title{
The Phylogenetic Regression of Grafen (1989)
}
\description{
The Phylogenetic Regression provides general linear model facilities for cross-species analyses, including hypothesis testing and parameter estimation. It uses branch lengths to account for recognised phylogeny (which makes the errors of more closely related species more similar), and the single contrast approach to account for unrecognised phylogeny (that a polytomy usually represents ignorance about which exact binary tree is true, and so one higher node should contribute only one degree of freedom to the test). One dimension of flexibility in the branch lengths is fitted automatically. Version 1.0.1 is a purely technical amendment of Version 1.0 to satisfy changed requirements for R packages.}
\details{
\tabular{ll}{
Package: \tab phyreg\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2018-04-05\cr
License: \tab GPL-2 | GPL-3\cr
}
You need a dataset of species data, and a phylogeny for those species that is in either (i) a series of taxonomic vectors, (ii) a "phylo" object, (iii) in newick format, or (iv) a single vector of the kind used internally in this package. For branch lengths, you can use the default "Figure 2" method of Grafen (1989), or specify the height of each node, or give a height to each level of the taxonomic vectors. Then you can test for H0 of one model against an HA of another, where the difference can be one or more x-variables, or interactions, or both, where the x-variables can be continuous or categorical. You may choose which species to include and exclude in each analysis. Missing data is handled automatically: you needn't do anything special, though a species is simply omitted from an analysis for which it has missing data in the response or in the control or test variables. It is assumed that missing data has the standard R value of NA.

In all linear models, the logic of each test involves controlling for some model terms while adding some test term(s). In simple cases such as ordinary multiple regression, the same analysis will give lots of tests and the user need only work out what the (often implicit) control and test terms are for each test. With the phylogenetic regression, only one test can be performed for a given analysis, and it is necessary to be explicit each time about the control terms and the test terms. This is because the single contrast taken across the daughters of one higher node depends upon the residuals in the control model.

Full details and examples are given under \code{\link{phyreg}}
}
\author{
	
Alan Grafen, with portions copied as follows.

(1) \code{read.newick} (used internally only) copied from Liam Rewell (see \code{\link{phyfromnewick}} for a more detailed acknowledgment)

(2) the definition of \code{ginv} has been copied from MASS (package comes with current R downloads; book is W.N. Venables and B.D. Ripley (2002) Modern Applied Statistics in S. (Fourth Edition), Springer -- \url{http://www.stats.ox.ac.uk/pub/MASS4}). This avoids having to \code{require} the whole package, though it may mean I have amend it if R and MASS make a simultaneous change in the low-level routines they use. \code{ginv} was copied and pasted from R 3.0.2 on 64-bit MacOS on 24th January 2014. It finds the generalised inverse of a matrix.

(3) code for dealing with model formulae (internal functions \code{merge.formulae.ag} and \code{merge.formulae.test.ag}) was adapted from code of Steven Carlisle Walker obtained from 
\url{https://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/} in January 2014.

Maintainer: Alan Grafen <alan.grafen@sjc.ox.ac.uk>
}
\references{
 Grafen, A. 1989. The phylogenetic regression. Philosophical Transactions of the Royal Society B, 326, 119-157.
 Available online at \url{http://users.ox.ac.uk/~grafen/cv/phyreg.pdf}. Some further information including GLIM and SAS implementations is available at \url{http://users.ox.ac.uk/~grafen/phylo/index.html}.
 }
\seealso{
	\code{\link{phyreg}}
}
\keyword{ package }
