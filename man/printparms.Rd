\name{Seeing the parameter values}
\alias{printparms}
\title{
Prints the current parameter values for \code{phyreg}
}
\description{
This is useful because \code{phyreg} has a mechanism for remembering your parameter values, so you don't have to specify them unless you're changing them, on each successive call. But this means you might forget, or want to check up after some puzzling output, just exactly what did I ask \code{phyreg}?}
\usage{
printparms()
}
%- maybe also 'usage' for other objects documented here.
\details{
Prints the current, remembered, parameters. It is in a bald print format, as this makes maintenance of the command trivial. It may also be useful to know that NULL and "NULL" as values for a parameter are treated by the program as the same.

Calling \code{printparms()} before your first call of \code{phyreg} should produce a warning message and no output, as only on the R session's first (successful) call of \code{phyreg} is a set of parameters remembered.
}
\author{
Alan Grafen}
\seealso{
\code{\link{phyreg}} tells you what the parameters mean, and has an example.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Parameter values }
