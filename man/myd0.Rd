\name{Example datasets}
\alias{myd0}
\alias{myd1}
\alias{myd2}
\alias{myd3}
\alias{extax}
\alias{newickstr}
\alias{exphylo}
\docType{data}
\title{
Example data
}
\description{
These datasets are supplied to allow examples to run in the \code{\link{phyreg}} help page. They all correspond to a hypothetical set of 100 species, with the same phylogeny and branch lengths.}
\usage{data(myd0)
data(myd1)
data(myd2)
data(myd3)
data(extax)
data(newickstr)
data(exphylo)}
\format{
  \code{myd0} to \code{myd3} are datasets with 100 observations on the following 7 variables.
  \describe{
    \item{\code{error}}{drawn from a standard Normal distribution for each species, but appropriately covarying across species.}
    \item{\code{X1-X3}}{continuous variables, each individually standard Normal, but appropriately covarying across species.}
    \item{\code{A, B}}{factors with four levels \code{1} to \code{4}}
    \item{\code{y}}{the response variable.}
  }
  \code{extax} is a data frame containing three taxonomic variables describing a phylogeny for 100 species
  \code{newickstr} is a character variable containing the same phylogeny as represented by \code{extax}, but in newick format. \code{exphylo} is a \code{phylo} object containing the same phylogeny as represented by \code{extax}. See Description for explanation of the heights. Note that a singleton higher node (only one daughter node) present in \code{newickstr} is retained in \code{exphylo}.
}
\details{
The formula calculating y from the independent variables was the same in each dataset, but the independent variables and the error were resampled for each dataset. The mean square error increases from \code{myd0} to \code{myd3}. By including "error" in an \code{lm()} (never possible in life, of course) you can obtain the coefficients and root mean square error used to construct the data. (\code{phyreg} will fail with no residual error, but without error the estimates are the same.) The categorical characters are evolved on the tree by taking a transition matrix (\code{tm}) to represent the transition probabilities after one unit of time, taking the matrix logarithm (\code{ltm<-logm(tm)}), and using the transition matrix \code{expm(bl * ltm)} for a branch length of duration \code{bl}. (Note that \code{expm()} and especially \code{logm()} are available in the package \code{expm}.)

For examples of use, see \code{\link{phyreg}}.
}
\seealso{The data is useful in calls to \code{\link{phyreg}}, with arguments specified as \code{data=myd0} etc and \code{taxmatrix=extax}. \code{newickstr} and \code{exphylo} are useful in calls to \code{phyfromnewick} and \code{phyfromphylo}.
}
\examples{
data(myd0, myd1, myd2, myd3, extax, newickstr, exphylo)
}
\keyword{datasets}
\keyword{phylogenetic data}

