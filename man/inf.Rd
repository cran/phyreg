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
inf(inpr, ..., oppf = 5, opdf = 1, oprho = 0, dfwarning = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inpr}{
An object of class \code{phyreglm} i.e. that has been produced from an assignment of the form \code{m<-phyreg(...)}.
}
  \item{\dots}{
To see the list of the possible parameters, call \code{inf()} itself with no parameters. They are arguments of \code{phyreg} that influence printing and storage of information. Four are shown below and can vary the way information is presented compared to the original call. See Details for information on the others.
}
  \item{oppf}{
If non-zero, the p-value and F-statistic are printed out at each call of \code{phyreg}, with \code{ndigits=oppf}. If zero, they are not printed out. \code{(inf)} will print these out for you if you include \code{"oppf"}.
}
  \item{opdf}{
If non-zero, a breakdown of degrees of freedom is printed out at each call of \code{phyreg}. \code{(inf)} will print this out for you if you include \code{"opdf"}.
}
  \item{oprho}{
If nonzero, details of rho will be printed with \code{ndigits=oprho}. If rho was fitted, the value of rho and the log-likelihood will be given, and how the search ended (lower boundary, or an internal maximum). \code{(inf)} will print this out for you if you include \code{"oprho"}.
}
  \item{dfwarning}{
If \code{TRUE} or \code{1}, the loss of degrees of freedom for phylogenetic reasons is detailed. \code{(inf)} will print this out for you if you include \code{"dfwarning"}.
}
}
\details{The \code{...} are used rather irregularly to allow values to be specified as character constants. These are the 13 final arguments of \code{phyreg()}, which specify material to be printed or to be stored. The printed material will be reprinted by \code{inf()}, saving the need to repeat the analysis. For saved material, including the character constant will instruct \code{inf()} to tell you whether that material was stored and both how to access it if it was stored and , if not, how to store it next time you use \code{(phyreg)}. You may include as many as these character constants as you like in one call of \code{inf()}. The "printing" constants are \code{"parmx"}, \code{"parmxz"}, \code{"opfunccall"}, \code{"means"}, and \code{"addDF"}, while the "saving" constants are \code{"linputs"}, \code{"sinputs"}, \code{"lmshortx"}, \code{"lmshortxz"}, \code{"lmlongx"}, \code{"lmlongxz"}, \code{"hinput"} and \code{"paper"}. To find what information is stored or printed by each character constant, see the help file for \code{phyreg} and look for the argument with the same name.
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
