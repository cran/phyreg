\name{Example datasets}
\alias{myd0}
\alias{myd1}
\alias{myd2}
\alias{myd3}
\alias{extax}
\alias{newickstr}
\docType{data}
\title{
Example data
}
\description{
Datasets \code{myd0} to \code{myd3} are data frames that contain a y-variable, three continuous x-variables \code{X1} to \code{X3}, and two factors \code{A} and \code{B}. They were simulated by a phylogenetic method, and also contain the error. The phylogeny was created from the taxonomic variables in \code{extax} -- the heights of the levels were assumed to be 0 for species, 1,2,3 for the three taxonomic variables, and 4 for the root.}
\usage{data(myd0)
data(myd1)
data(myd2)
data(myd3)
data(extax)
data(newickstr)}
\format{
  \code{myd0} to \code{myd3} are datasets with 100 observations on the following 7 variables.
  \describe{
    \item{\code{error}}{drawn from a standard Normal distribution for each species, but of course correlated appropriately for the phylogeny}
    \item{\code{X1-X3}}{continuous variables}
    \item{\code{A, B}}{factors with four levels \code{1} to \code{4}}
    \item{\code{y}}{the response variable.}
  }
  \code{extax} is a data frame containing three taxonomic variables describing a phylogeny for 100 species
  \code{newickstr} is a character variable containing the same phylogeny as represented by \code{extax}, but in newick format -- see Description for explanation of the heights.
}
\details{
The formula calculating y from the independent variables was the same in each dataset, but the independent variables and the error were resampled for each dataset. The mean square error increases from \code{myd0} to \code{myd3}. By including "error" in the model (never possible in life, of course) you can obtain the coefficients and root mean square error used to construct the data. The categorical characters are evolved on the tree by taking a transition matrix (\code{tm}) to represent the transition probabilities after one unit of time, taking the matrix logarithm (\code{ltm<-logm(tm)}), and using the transition matrix \code{expm(bl * ltm)} for a branch length of duration \code{bl}. (Note that \code{expm()} and especially \code{logm()} are available in the package \code{expm}.)

For examples of use, see \code{\link{phyreg}}.
}
\seealso{The data is useful in calls to \code{\link{phyreg}}, with arguments specified as \code{data=myd0} etc and \code{taxmatrix=extax}. \code{newickstr} is useful in calls to \code{phyfromnewick}.
}
\examples{
data(myd0, myd1, myd2, myd3, extax, newickstr)
}
\keyword{datasets}
\keyword{phylogenetic data}

