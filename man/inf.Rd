\name{inf}
\alias{inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provides information about the stored output of a call to phyreg()
}
\description{
After an instruction \code{m<-phyreg(...)}, you can type \code{inf(m,...)} to provide information about \code{m}. You can see output again, or that you didn't ask for from the initial call of \code{phyreg()}. And you can find out whether optional storage took place of the different possible outputs, if not, how to ask for it, and if so, how to access it.
}
\usage{
inf(inpr, ..., oppf = 5, opdf = 1, oprho = 0, dfwarning = 0, plot=0, others)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inpr}{
An object of class \code{phyreglm} i.e. that has been produced from an assignment of the form \code{m<-phyreg(...)}.
}
  \item{\dots}{
To see the list of the possible parameters, call \code{inf()} itself with no parameters. They are arguments of \code{phyreg} that influence printing and storage of information. They should be supplied in quotation marks e.g. \code{inf(m,"sinputs","means")}.
}
  \item{oppf}{
If non-zero, see the p-value an F-ratio for the test conducted in the \code{phyreg()} call, with \code{ndigits=oppf}.
}
  \item{opdf}{
If non-zero, see a breakdown of the degrees of freedom in the test conducted by the \code{phyreg()} call.
}
  \item{oprho}{
If non-zero, see rho-related information from the \code{phyreg()} call, with \code{ndigits=oprho}.
}
  \item{dfwarning}{
If non-zero, see warnings about the loss of degrees of freedom for phylogenetic reasons in the test conducted by the \code{phyreg()} call. See Details under \code{\link{phyreg}} for an explanation.
}
  \item{plot}{
If non-zero, this produces the plot that \code{phyreg} can produce at the time of the call, so see \code{\link{phyreg}} for an explanation of the plot itself.}
  \item{others}{A list of further values are permitted, all need to be in quotation marks. There is no need to specify them as, for example, \code{others="means"}, and they can be given directly as separate arguments to \code{inf()}. They are all parameters of \code{phyreg}. First, \code{"oppf"},\code{"opdf"},\code{"oprho"}, \code{"dfwarning"}, \code{"plot"} can be specified in this way, in which case they receive default values of 5, 1, 5, 1, and 1, respectively. A second batch control the printing out of values, namely \code{"parmx"}, \code{"parmxz"},\code{"opfunccall"},\code{"means"}, and \code{inf} will print these for you without having to rerun \code{phyreg}. Finally, the storage of information is determined by  \code{"lmshortx"}, \code{"lmshortxz"}, \code{"lmlongx"}, \code{"lmlongxz"}, \code{"paper"}, \code{"linputs"}, \code{"sinputs"}, \code{"hinput"}. \code{inf} will tell you whether the listed items were stored in that particular call, how to store it next time if it wasn't this time, and how to access it if was stored this time.}
}
\details{
Examples are given under \code{\link{phyreg}}.
}

\author{
Alan Grafen
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{phyreg}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ information }
\keyword{ phyreg }% __ONLY ONE__ keyword per line
