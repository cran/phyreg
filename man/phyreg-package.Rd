\name{phyreg-package}
\alias{phyreg-package}
\docType{package}
\title{
The Phylogenetic Regression of Grafen (1989)
}
\description{
The Phylogenetic Regression provides general linear model facilities for cross-species analyses, including hypothesis testing and parameter estimation. It uses branch lengths to account for recognised phylogeny (which makes the errors of more closely related species more similar), and the single contrast approach to account for unrecognised phylogeny (that a polytomy usually represents ignorance about which exact binary tree is true, and so one higher node should contribute only one degree of freedom to the test). One dimension of flexibility in the branch lengths is fitted automatically. It incorporates "independent contrasts" and "phylogenetic generalised least squares" of current theory -- see Details below for a discussion.}
\details{
\tabular{ll}{
Package: \tab phyreg\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2014-02-08\cr
License: \tab GPL-2 | GPL-3\cr
}
This package fits into my understanding of current theories as follows. Grafen (1989) begins with a Phylogenetic Generalised Least Squares model (which he calls "the standard regression", as it depends on standard statistical theory), and converts it into independent contrasts, proving in Theorem 1 they are completely equivalent (called the "long regression", with one datapoint for each branch of the tree). However, this regression is \emph{too} long! He adopts the Ridley principle (Ridley 1983) that each higher node should contribute just one datapoint, and to do this finds a single contrast across the daughters of a higher node. It is sometimes asserted that this choice of a single contrast, like other suggested choices, is arbitrary (for example, Garland and Diaz-Uriarte 1999) -- however, sections 3(b) and 3(c) of Grafen (1989) provides lines of argument to show that the choice there is not only natural, but the right one. I am not aware of any engagement with the substance of that argument, and so to my mind it remains entirely unrefuted. The set of single contrasts forms the "short regression". Theorem 2 of Grafen (1989) provides a formal mathematical statement that the F-ratio in the short regression actually does have an F-distribution with the relevant degrees of freedom. Thus, that test is much more precisely defined and defended than the suggestion of Purvis and Garland (1993) of taking the F-ratio as calculated from the independent contrasts, but comparing it to critical F-ratios with reduced denominator degrees of freedom. Methodology often involves fudges where they are necessary -- my claim is that, contrary to what appears to be the accepted view in the literature, no fudge at that level is needed here.

This package provides all the basic output for all the regressions we've discussed. As the PGLS and independent contrasts regressions are equivalent, you can find the residual sums of squares, parameter estimates, t-tests for them, etc, in the output from the "long regression". The single contrasts regression is the test provided by the package, and its detailed output is under the "short regression". One element of the PGLS output is not provided by the long regression, namely the fitted values, so the output of \code{phyreg} contains them as \code{pglsFVx} and \code{pglsFVxz}, for the control and control+test models, respectively. From them you can calculate the residuals as y-variable minus fitted values. It is important, however, to say that the only hypothesis testing that is justified is the F-ratio provided by \code{phyreg} -- the others you can obtain or construct are invalid, and the only use I can see for them is to show quite how wrong they are.

If it sounds like I'm saying "I was right 25 years ago, and am still right today", you have understood me exactly, but this is not a matter of mere assertion. The key question is whether the single contrasts I extract are the right ones. I am not at the time of writing aware of any substantial discussion of that point -- perhaps it was too technical at the time. You may like to read sections 3(b) and 3(c) of Grafen (1989), and also Grafen (1992), which explicitly tackles the question of whether other contrasts might work, and adds some further points. For example, if the test is for adding two independent variables X and Z, then we may well desire a test that would give the same answer if we instead proposed to add X-Z and X+Z, as they contain the same information -- contrasts based on test variables cannot do this, while it is automatic under the phylogenetic regression. Unless all those arguments are overturned, to my mind the phylogenetic regression is indeed the (only) right way to do linear models on comparative data. Fortunately, it is now available in R!

There is one unusual feature of the method that I draw attention to here. In all linear models, the logic of each test involves controlling for some model terms while adding some test term(s). In simple cases such as ordinary multiple regression, the same analysis will give lots of tests and the user need only work out what the (often implicit) control and test terms are for each test. For example, in usual regressions, Type I sums of squares test for each variable, controlling for all previous variables, while Type III sums of squares test for each variable, controlling for all the others. With the phylogenetic regression, only one test is performed for a given analysis, and it is necessary to be explicit each time about the control terms and the test terms. This is because the single contrast taken across the daughters of one higher node depends upon the residuals in the control model, and so one analysis must always have the same control variables in tests it provides. In principle, then, we could have a single fixed set of control variables, but a number of different sets of test variables -- but currently, only one set of test variables is handled. (A major gain in efficiency can be obtained for this situation \emph{on the assumption you have no missing values}. Obtain rho for an analysis with a given set of control variables, and then set rho in the arguments of subsequent calls to phyreg with the same control variables. This works because rho is fitted just to the control model, and most of the execution time is spent fitting rho. But the fitted value of rho will also depend on which species are included, hence the caveat about missing values.)

One key conceptual difference between Felsenstein's 1985 paper and my 1989 paper is that Felsenstein assumed that the traits themselves undergo Brownian motion over evolutionary time, whereas I did my analysis on the basis that only the error term does so. (This point is rediscovered every now and then.) It is important because (i) it is a much smaller assumption (ii) it allows independent variables to be categorical (iii) it means that as we add independent variables to a model, the correct branch lengths may well change, because the error before you add X includes X, while afterwards it doesn't, so removing a variable that explains a lot of variation between major groups may well mean we should be altering our branch lengths to increase the relative lengths of those nearer the species end. Thus there is no hard and fast principle about branch lengths reflecting durations between splits, or indeed anything else - they must treated pragmatically. I also don't see any merit in alternative evolutionary processes (such as Ornstein-Uhlenbeck) when applied to the error -- however appropriate they may or may not be for the traits themselves.

Full details and examples are given under \code{\link{phyreg}}
}
\author{
	
Alan Grafen

Maintainer: Alan Grafen <alan.grafen@sjc.ox.ac.uk>
}
\references{
 Felsenstein, J. 1985. Phylogenies and the comparative method. American Naturalist 125:1-15.
 
 Garland, T., Jr, and R. Diaz-Uriarte, 1999. Polytomies and phylogenetically independent contrasts: an examination of the bounded degrees of freedom approach. Systematic Biology 48: 547-558.
 
 Grafen, A. 1989. The phylogenetic regression. Philosophical Transactions of the Royal Society B, 326, 119-157. Available online at \url{http://users.ox.ac.uk/~grafen/cv/phyreg.pdf}. Some further information, including GLIM (Trademark -- see http://www.nag.com) and SAS (Trademark -- see http://www.sas.com) implementations, is available at \url{http://users.ox.ac.uk/~grafen/phylo/index.html}.
 
 Grafen,. A. 1992. The uniqueness of the Phylogenetic Regression. Journal of theoretial Biology 156: 405-423.
 
 Purvis, A. and T. Garland, Jr. 1993. Polytomies in comparative analyses of continuous data. Systematic Biology 42: 569-575.
 
 Ridley, M. 1983. The explanation of organic diversity. Oxford: Clarendon Press.}
\seealso{
	\code{\link{phyreg}}
}
\keyword{ package }
