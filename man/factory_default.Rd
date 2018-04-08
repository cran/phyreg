\name{factory_default}
\alias{factory_default}
\alias{new_default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default settings for "minor" parameters
}
\description{
There are many minor parameters of \code{phyreg} that affect output, storage of statistical output, and less important aspects of model fitting. Major and minor parameters are remembered within a session, but the minor parameters are set by default when the package is loaded, and can be returned to by including \code{reset=TRUE} as an argument of \code{phyreg()}. You can choose to store the current remembered minor parameter values as those defaults by calling \code{new_default()}, and you can choose to return the default values to those distributed with the package by calling \code{factory_default()}.
}
\usage{
factory_default()
new_default()
}
%- maybe also 'usage' for other objects documented here.
\details{
The unstored, major, parameters are the \code{control} and \code{test} parameters, and the data parameters, namely \code{data}, \code{phydata}, \code{taxmatrix} and \code{heightsdata}. All the parameters are listed under \code{\link{phyreg}}, and examples of use are also given there.
}
\author{
Alan Grafen
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{phyreg}}
}

\keyword{ default parameters }
