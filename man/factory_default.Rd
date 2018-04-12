\name{factory_default}
\alias{factory_default}
\alias{new_default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default settings for "minor" parameters
}
\description{
These functions change the effect of including \code{reset=TRUE} as an argument to \code{phyreg()}. After a call of \code{factory_default()}, the reset will be to the original values hard-wired in the package's code. If you call \code{new_default()}, the values of the minor parameters in place at that moment become the reset values for all uses of \code{reset=TRUE} until a future call of one of these functions.
}
\usage{
factory_default()
new_default()
}
%- maybe also 'usage' for other objects documented here.
\details{
The unstored, major, parameters are the \code{control} and \code{test} parameters, and the data parameters, namely \code{data}, \code{phydata}, \code{taxmatrix} and \code{heightsdata}. The remaining parameters are minor. All the parameters are listed under \code{\link{phyreg}}, and examples of use are also given there. Note that changing a parameter whether major or minor leaves the newly set value of that parameter as the default for future calls of \code{phyreg()} until either the parameter is changed again, or \code{reset=TRUE} is included as an argument of \code{phyreg()}.
}
\author{
Alan Grafen
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{phyreg}}
}

\keyword{ default parameters }
