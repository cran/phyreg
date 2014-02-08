\name{phyreg}
\alias{phyreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs a phylogenetic regression 
}
\description{
This function performs a single test in a linear models framework of an alternative hypothesis against a null hypothesis in cross-species data, taking account of the phylogeny of the species, implementing the phylogenetic regression of Grafen (1989).

The hypotheses are specified in R format as model formulae or model terms. The phylogeny can be specified either using an internal "phy" format or a matrix of taxonomic variables -- conversion functions are supplied from the standard newick and phylo formats. Node heights can be specified to determine the phylogenetic covariance structure of the error, or left to default values. The user can specify a value of rho (a parameter that determines whether branch lengths are bigger near the root or near the tips), or ask that it be fitted by maximum likelihood.

To get acquainted, bother only with the first three arguments (control model, test terms, and dataset) and then either the fourth or fifth to specify the phylogeny in one or other of two ways. The rest can wait!

The default outputs provides a p-value and F-ratio, and repeats the control and control+test models. Much more information is available, though much of that needs to be interpreted with care.
}
\usage{
phyreg(control, test, data, subset, phydata, taxmatrix, heightsdata, 
rho = -1, lorho = 0.3, hirho = 0.6, errrho = 0.02, minrho = 1e-04, 
tolerance = 1e-06, oppf = 5, opdf = 0, parmx = 0, parmxz = 0, 
opfunccall = 0, addDF = 0, linputs = FALSE, sinputs = FALSE, means = FALSE, 
lmshortx = FALSE, lmshortxz = FALSE, lmlongx = FALSE, lmlongxz = FALSE, 
hinput = FALSE, paper = FALSE, dfwarning = TRUE, oprho = FALSE, plot=FALSE, 
reset=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{Arguments are remembered from one call to the next -- so you need to specify only the ones that change. You can disable this by including \code{reset=TRUE} in the your call of \code{phyreg}.
  \item{control}{
The null hypothesis model, which should be a formula with a response variable. If there are no variables to control for, specify \code{y~1}. (\code{y~0} is not permitted.)
}
  \item{test}{The test terms. To control for \code{y~b}, and test for \code{a*x}, you can specify (i) \code{test="a*x"} or (ii) \code{test=.~b+a*x} (Note the LHS *must* be a dot for this "update-style" specification), or (iii) \code{test=~+a*x}. (Note the "+" \emph{and} the absence of the response variable in this final, "adding terms", form.)
}
  \item{data}{A mandatory dataframe containing the variables in the model formulae}
  \item{subset}{
a vector of the same length as the data, containing 0/1 for exclusion or inclusion of a species. If unspecified, the internal default is to include all species. To return to being unspecified, say \code{subset=NULL}
}
  \item{phydata}{
a vector with the phylogeny in the internal package format (containing for each non-root node the ID number of its parent). Either \code{phydata} or \code{taxmatrix} must be supplied. To unspecify, say \code{phydata=NULL}. A phylogeny in newick or phylo format can be converted into the internal format, preserving height information, with \code{phyfromnewick} or \code{phyfromphylo}.
}
  \item{taxmatrix}{
a dataframe with a set of taxonomic vectors to specify the working phylogeny. There should be no species column, and the others should start with low levels (e.g. genus) and end with the highest (perhaps order). If \code{taxmatrix} is specified, you can alternatively specify \code{heightsdata} with a vector containing one height for each level you've given, plus a height for the root, or if that is \code{NULL} the default "Figure 2" heights of Grafen (1989) will be used. Either \code{taxmatrix} or \code{phydata} (-- not both) must be supplied. Despecify with \code{taxmatrix=NULL}.
}
  \item{heightsdata}{
If not specified, the height of nodes will be calculated using the default "Figure 2" method of Grafen (1989). If a height is specified for each node, then the vector will be one element longer than \code{phydata}, as the root has a height but no parent. If \code{taxmatrix} is specified, then \code{heightsdata} can also be specified as a vector of length one plus the number of taxonomic vectors, in which case each node is given a height corresponding to its level in the taxonomy -- the extra element is the height of the root, and the species are assumed to have a height of zero. A parent node must be higher than its daughter -- zero branch lengths are not allowed, but are equivalent to combining the daughters of the two nodes under a single node.
}
  \item{rho}{
If zero or negative, \code{rho} will be fitted by maximum likelihood. If positive, that value will taken for rho and no fitting will be done. Once the node heights are scaled to lie between 0 and 1, each height is raised to the power rho before being used to determine the error structure in the model. See Grafen (1989) for details.
}
  \item{lorho}{
The starting lower bound for the search for rho
}
  \item{hirho}{
The starting upper bound for the search for rho
}
  \item{errrho}{
The search for rho will stop when it is known to lie in an interval of length \code{errrho}
}
  \item{minrho}{
The search for rho will stop when the whole current search interval lies below \code{minrho}. A zero value of rho is problematic for the program.
}
  \item{tolerance}{
This tolerance decides whether there is enough variability among the daughters of a higher node to justify putting it into the short regression. This is because the residuals in the long regression on the control variables are used to calculate a contrast; if they are very small, then the contrast is effectively being determined by machine rounding errors, and if they are zero, then all contrasts should also be zero and so the higher node will be irrelevant in a regression through the origin.
}
  \item{oppf}{
If non-zero, the p-value and F-statistic are printed out at each call of \code{phyreg}, with \code{ndigits=oppf}. If zero, they are not printed out.
}
  \item{opdf}{
If non-zero, a breakdown of degrees of freedom is printed out at each call of \code{phyreg}.
}
  \item{parmx}{
If non-zero, the parameters in the long regression with the control formula are printed out at each call of \code{phyreg}.
}
  \item{parmxz}{
If non-zero, the parameters in the long regression with the control+test formula are printed out at each call of \code{phyreg}. Note these should be regarded with caution, as the value of rho was fitted just on the control formula. For better estimates, run \code{phyreg} again with the old control+test as the new control, and set \code{opparmx=1}. But \code{opparmxz} will usually give a good rough guide.
}
\item{plot}{A plot is produced based on the output of the short regression, which has one point for each higher node included. The data for this plot are contained after \code{x<-phyreg(...)} as \code{x$details$influence$RCT} and \code{x$details$influence$DR}. The horizontal axis, \code{RCT}, is the residual of the point after fitting the control+test model. The vertical axis, \code{DR}, contains the improvement in that point through including the test terms, measured as the difference between the squares of the residuals. In each case, the residuals are adjusted for degrees of freedom, i.e. divided by the residual degree of freedom of the control or control+test model as appropriate. The variables are then scaled to make them easier to read and interpret. \code{DR} is scaled so it adds to 100, and each point then indicates the percentage of the improvement in SS that is due to a higher node (possibly negative if it fits less well after the test variables are added). \code{RCT} is scaled so the mean square error is \code{1}, and so Normal errors would give 95\% of points within 1.96 of the mean (but remember the regression is through the origin so the average of \code{RCT} does not have to be zero). The significance of the plot is that extreme points on the horizontal axis fit badly after adding the test variables, while those near zero fit well. Points near zero on the vertical axis have been little changed in their contribution to badness-of-fit by adding the test variables, while highly positive points have had a lot of their deviation explained by the test variables.}
  \item{opfunccall}{
If non-zero, the function call is printed out
}
  \item{addDF}{
Advised to leave well alone. Included for backwards compatibility with the original GLIM (Trademark) version. It was originally but misguidedly intended to make a correction to the total degrees of freedom for the test in recognition of the fitting of rho. However, the fitting of rho is now regarded as affecting the numerator and denominator SS equally, and so not to require the correction.
}
  \item{linputs}{
If \code{TRUE} or \code{1}, the inputs to the long regression will be stored. After \code{x<-phyreg(...,linputs=1)}, it will be found that \code{x$linputs$y}, \code{x$linputs$designxz} and \code{x$linputs$w} will contain the response, design matrix, and weights for the long regression, respectively. Use with care. The whole point of the the problem of \emph{unrecognised} phylogeny is that the standard errors from this regression cannot be trusted, though with a dependable binary phylogeny, the output can be accepted at face value.
}
  \item{sinputs}{
If \code{TRUE} or \code{1}, the inputs to the short regression will be stored. After \code{x<-phyreg(...,sinputs=1)}, it will be found that \code{x$sinputs$y}, \code{x$sinputs$designxz} and \code{x$sinputs$w} will contain the response, design matrix, and weights for the short regression, respectively. Use with care. Not only can the standard errors from this regression not be trusted, but the parameter estimates are biassed!
}
  \item{means}{
If \code{TRUE} or \code{1}, the phylogenetically-weighted means will be printed for the response and for each variable in the design matrix of the control+test model (in the original species dataset, which are the same as the weighted means in the long regression).
}
  \item{lmshortx}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,lmshortx=1)}, the \code{lm.wfit()} output from the short regression with the control formula will be stored in \code{x$lmshortx}. Use with care -- the parameter estimates are biassed...
}
  \item{lmshortxz}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,lmshortxz=1)}, the \code{lm.wfit()} output from the short regression with the control+test formula will be stored in \code{x$lmshortxz}. Use with care -- the parameter estimates are biassed...
}
  \item{lmlongx}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,lmlongx=1)}, the \code{lm.wfit()} output from the long regression with the control formula will be stored in \code{x$lmlongx}. The SEs cannot be trusted for non-binary phylogenies.
}
  \item{lmlongxz}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,lmlongxz=1)}, the \code{lm.wfit()} output from the long regression with the control+test formula will be stored in \code{x$lmlongxz}. The SEs cannot be trusted for non-binary phylogenies.
}
  \item{hinput}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,hinput=1)}, the heights used will be stored in \code{x$hinput}. If you use \code{taxmatrix=...}, this is a way to obtain the internal "phy" version of the heights for each higher node. These heights have \emph{not} been modulated by rho i.e. they are the starting heights before rho is applied.
}
  \item{paper}{
If \code{TRUE} or \code{1}, then after \code{x<-phyreg(...,paper=1)}, various matrices that appear in the appendix of the source paper (Grafen 1989) will be stored under \code{x$paper}. Type \code{names(x$paper)} to see which ones!
}
  \item{dfwarning}{
If \code{TRUE} or \code{1}, the loss of degrees of freedom for phylogenetic reasons is detailed. It is a "warning", so if none are lost, nothing is printed, and it doesn't occupy space in your output.}
  \item{oprho}{
If nonzero, details of rho will be printed with \code{ndigits=oprho}. If rho was fitted, the value of rho and the log-likelihood will be given, and how the search ended (lower boundary, or an internal maximum).
}
\item{reset}{If \code{TRUE} or {1}, the parameters are all reset to the start-of-session values before applying the other parameters you specify in the same call of \code{phyreg}. Useful if you want to start a new project, or just be sure there's no previous setting that is lurking and affecting your results.}
}
\details{One key point for using the program is that each node in the phylogeny has an ID number, with which it is identified. For a species, the number is given by its position in the data frame with the data itself, which is assumed to have species in the same order as in the taxonomic variables data frame, if that is specified, and the same order as \code{phydata} if you use that way of specifying a phylogeny. When converting from a newick-style phylogeny, again the order of the species in the newick format is carried over to the \code{$phy} variable, and so determines their ID numbers. If \code{phyfromphylo} is used, the code numbers of the species in the "phylo" object determine the order they are placed in the internally formatted "phy" vector, so the original code numbers become their ID numbers. These species ID numbers matter because they must be aligned with the dataset, and are fixed no matter whether that species or another is excluded by the \code{subset} parameter of \code{phyreg}, or for missing data. Higher nodes have ID numbers too. In every phylogeny, each node's parent has a higher ID number than it does, except the root, which has no parent. The numbers for higher nodes \emph{do} vary when species are omitted for any reason. A function is provided to identify higher nodes for you, by listing either their daughters or all the species below them -- see \code{\link{who}}.
	
The program is a linear model, and so allows the independent variables to be continuous or categorical, and allows models to contain interactions. R's conventions for model formulae are followed, so see \code{help(formula)} for the symbols for interaction and so on, and \code{help(factor)} for how to declare a variable as categorical (the default is continuous).

Missing data is handled automatically: you needn't do anything special, though a species is simply omitted from an analysis for which it has missing data in the response or in the control or test variables. It is assumed that missing data has the standard R value of \code{NA}.

You can choose which species to include or omit in a given analysis using the \code{subset} argument, set to 0 or 1 for omission or inclusion. Clear by asserting \code{subset=NULL}, and all species will be included, subject to missing data (see just above).
	
The main feature of phylogenetic data that the Phylogenetic Regression deals with differently from other methods is \emph{unrecognised phylogeny}. It operates on the principle that each higher node should provide one datapoint to a valid test, and not more, and does so by choosing just one linear contrast across the daughters of each higher node. (This is different from Purvis and Garland (1993)'s suggestion of calculating a test allowing more datapoints per higher node, and then looking up the critical value as though you hadn't.) The contrast coefficients come from the residuals of the long regression on the control model. For full mathematical and conceptual justifications of this choice, see Grafen (1989, 1992). With binary phylogenies, of course, there is no unrecognised phylogeny. A parameter of the branch lengths is fitted automatically, unless the user wishes to impose a value, which allows the strength of phylogenetic association to be make weaker or stronger. Simulation studies in Grafen (1989) show that the method has good properties, and also that ignoring unrecognised phylogeny can create serious problems.

For a discussion of the place of the Phylogenetic Regression in current theory and methods, see the \code{\link{phyreg-package}}.

On some occasions, degrees of freedom are lost "for phylogenetic reasons" (Grafen 1989, section 3(e)). A whole node may be lost to the final test if the residuals of its daughter nodes are all zero in the long regression. This can happen for various reasons, most often when (i) the response is in fact binary, and so there is no variation in it below a node, or (ii) a categorical variable has so many values restricted to one part of the tree that a subset of its parameter values can adjust to render all the residuals zero in that part of the tree. That is called a node being lost in the denominator. The other possibility is that once the contrasts have been taken across each higher node, the design matrix for the model has lower rank than it did before, which is called losing a degree of freedom in the numerator (it is transferred to the denominator). You can choose to be warned when degrees of freedom are lost for phylogenetic reasons (use \code{dfwarning=1}, or to see a whole breakdown of degrees of freedom including any lost (use \code{opdf=1}). If nodes are lost in the denominator, their ID numbers are stored in the output \code{$details$missingnodes}, though note this also contains an initial \code{0} for programming convenience. These phylogenetic degree of freedom issues need to be handled properly to provide a valid test.

The data variables (\code{data}, \code{phydata}, \code{heightsdata}, \code{taxmatrix}) are stored by name, and so if you change a dataset, and then call \code{phyreg}, those changes will be incorporated, even when that argument is just carried on by default from a previous call. The other variables, including \code{subset}, are stored by value, and so remain the same until set again in \code{phyreg}.
	
The data for long and short regressions, and the \code{lm.wfit()} output are provided on request, but must be used with great caution. The whole point of the phylogenetic regression is that neither type of regression provides results that can be taken at face value, except the long regression under binary phylogenies. Use at your own risk! See Grafen (1989) for details.

Some simulated data is included to facilitate the examples below -- see \code{\link{myd0}}
}
\value{
 \item{H0model}{The control model}
 \item{HAmodel}{The control+test model}
  \item{spu}{A vector with a 0/1 for each species to indicate whether it was used (1) or not (0), depending on the argument \code{subset} and on missing values}
 \item{nomspuse}{The number of species omitted because of missing values (not counting those already omitted because of the argument \code{subset})}
 \item{longrss}{The residual sum of squares in the long regression with the control formula }
 \item{missingnodes}{A vector of the ID numbers of nodes omitted for phyogenetic reasons (see Details above). An extra zero is included for programming convenience.}
 \item{shornode}{The ID numbers of the higher nodes in the supplied phylogeny that are used in the short regression, so \code{shornode[5]} gives the ID number of the node for the 5th datapoint.}
 \item{details}{contains various numbers needed for output, namely \code{$poff} (the p-value), \code{$testf} (the F-ratio), \code{$testnumdf} (the numerator DF), \code{$testdendf} (the denominator DF), \code{$nspec} (total number of species, included and excluded, so it equals \code{nrow(input_dataframe))}, \code{$nommiss} (the number of species omitted because of missing values (not counting those already omitted because of the argument \code{subset})), \code{$nspecactive} (number of species included in the analysis), \code{$nomspuse} (number of species omitted because of the argument \code{subset}), \code{$nphyhinodes} (the number of non-species nodes in the supplied phylogeny (ie. with all species included)), \code{nhinodes} (number of higher nodes in the phylogeny as used, which may be reduced by missing data or by \code{subset}), \code{$hilostomsp} (the number of higher nodes lost through species omitted because of the argument \code{subset}, \code{$hilostvar} (number of higher nodes lost to the short regression through lack of variability in the long regression -- see "Phylogenetic Degrees of Freedom" in Details), \code{$shortotdf} (total DF in short regression), \code{$dflx} (rank of the long regression with the control model), \code{$dfxlost} (loss of numerator DF in short regression -- see "Phylogenetic Degrees of Freedom" in Details), \code{$dfsx} (rank of short regression with control model), \code{$testlongdf} (\code{=$dflxz-$dflx}, the degrees of freedom for the test terms in the long regression), \code{$testlost} (degrees of freedom for the test variables, lost for phylogenetic reasons -- see "Phylogenetic Degrees of Freedom" in Details)), \code{$shortcondf} (control degrees of freedom in the short regression), \code{$addDF} (retains the argument \code{addDF} as it's needed in calculations), \code{$rho} (the value of rho used in the final test), \code{$lik} (the log-likelihood of that value of rho), \code{$edge} (coded -1 for error, 2 for rho set by the user, and, with a fitted rho, 0 for an internal maximum of the likelihood, and 1 for a maximum at the lower edge of the search region (specified by \code{minrho})), \code{$missingnodes} (the numbers of the higher nodes omitted for lack of variability in the short regression, but an extra 0 at the start -- see Details), and finally a data frame \code{$influence} containing three variables about the addition of the test variables to the short regression (\code{NN} identifies the higher node with its ID number, and \code{$RCT} and \code{$DR} are measures for each node explained under the \code{plot} argument above -- use the function \code{\link{who}} for help pinning down what each higher node is).}
 \item{means}{The phylogenetically weighted means of the response variable and the design matrix for the control+test model}
 \item{parmx}{Parameter estimates from the long regression with the control model}
 \item{parmxz}{Parameter estimates from the long regression with the control+test model}
 \item{funccall}{The function call with which you invoked \code{phyreg}}
 \item{fullphy}{The full phylogeny for all species. This will just be \code{phydata} if you specified the phylogeny that way, but otherwise has been converted from taxonomic vectors.}
 \item{usedphy}{The phylogeny for the species included in the analysis. There is an entry for every species, and each included species has its original ID, and every omitted species has zero. There may well be fewer higher nodes in \code{usedphy} than in \code{fullphy}, because only higher nodes with two or more daughters are retained. This vector is useful if you want to study the matrices from the paper (see the argument and value \code{paper}) as many are indexed by \code{usedphy}. All node numbers reported by the program elsewhere are the original IDs, and so \code{usedphy} is important only for internal and/or technical reasons.}
 \item{originaIDs}{This tells you, for each node in \code{usedphy}, which node in \code{fullphy} it corresponds to. All the usual input already uses those \code{fullphy} names, so you would need this only if you were going exploring in \code{lmshortx(z)} or studying \code{usedphy}.}
 \item{sinputs}{(optional) the inputs to the short regression namely \code{sinputs$y}, \code{sinputs$design} and \code{sinputs$w} for the response, design matrix and weights. Choose to store using the argument of the same name.}
 \item{linputs}{(optional) the inputs to the long regression namely \code{linputs$y}, \code{linputs$design} and \code{linputs$w} for the response, design matrix and weights. Choose to store using the argument of the same name.}
 \item{lmlongx}{(optional) The \code{lm.wfit()} output from the long regression with the control model. Choose to store using the argument of the same name.}
 \item{lmlongxz}{(optional) The \code{lm.wfit()} output from the long regression with the control+test model. Choose to store using the argument of the same name.}
 \item{lmshortx}{(optional) The \code{lm.wfit()} output from the short regression with the control model. Choose to store using the argument of the same name.}
 \item{lmshortxz}{(optional) The \code{lm.wfit()} output from the short regression with the control+test model. Choose to store using the argument of the same name.}
 \item{pglsFVx}{These are fitted values for each species based on the control model only. They are calculated for all species with relevant data, and not just for those included in the model. The fitted values are those for (i) the phylogenetic generalised least squares model and (ii) the "standard regression of Grafen (1989)". They are related to (iii) the long regression of this package, as FVx = mean of y + sum over design variables of (coefficient for variable times (variable - its mean)). See \code{means} among the other values. You can calculate residuals as \code{residuals <- y - pglsGVx} }
 \item{pglsFVxz}{The same as \code{pglsFVx}, \emph{mutatis mutandis}, except the fitted values are based on the control+test model.}
}
\references{
	The source paper is Grafen, A. 1989. The phylogenetic regression. Philosophical Transactions of the Royal Society B, 326, 119-157. Some further information, including GLIM (Trademark -- see http://www.nag.com) and SAS (Trademark -- see http://www.sas.com) implementations, is available at \url{http://users.ox.ac.uk/~grafen/phylo/index.html}.
}
\author{
Alan Grafen
}

\seealso{
The package helpfile is at \code{\link{phyreg-package}}. The function \code{\link{inf}} displays information about stored \code{phyreg} output. \code{\link{phyfromnewick}} and \code{\link{phyfromphylo}} convert phylogenies from two different standard formats into the internal format for the package. For simulated data see \code{\link{myd0}}. A function \code{\link{printparms}} shows you the current remembered set of parameters for \code{phyreg}. To identify higher nodes from their ID numbers, use \code{\link{who}}.
}
\examples{
## First get some data
data(myd0,myd1, myd2, myd3, extax)
##
## Then do our first analysis
phyreg(y~X1,"A",myd0,taxmatrix=extax)
## and test instead for "B", noticing that only the changed parameter need be given
phyreg(test="B")
## and we do more complicated analysis involving an interaction with an existing term
phyreg(y~X1+X2,"A*X1")
##
## Now we choose to see the output relating to rho and to how the degrees of freedom 
##   are determined, and we also wish to see the means for each variable, and the  
##   parameters from the long regression on control+test variables
phyreg(oprho=6, opdf=1, means=1,parmxz=1)
## To illustrate inf, we store the results of an analysis in m1
m1<-phyreg(y~A,"X1+X2")
## Note we still get the extra output from the previous call, because those parameters 
##   too are remembered within a session. But we can see it again, whether or not we
##   saw it the first time, with inf. inf reminds you if you forget quite how to use it
inf()
inf(m1)
inf(m1,oppf=3)
inf(m1,oppf=7, oprho=5)
inf(m1, oppf=5, "means", "parmx")
inf(m1,"sinputs","shortx")
## The final call asks for things m1 doesn't have because it wasn't stored at the time. 
## inf tells you about remembered output. To find out about your own input you can
printparms()
##  This is the set of remembered internal parameters for phyreg.
##
##  The phylogeny has so far been determined by a data frame of taxonomic variables
##   in the argument taxmatrix. If we have the phylogeny available in newick style, we can  
##   convert to the interal format, and then use that instead. Fortunately, one is provided.
##   Note it is good form to unset the other method of specifying a phylogeny (which is 
##   being remembered by the package) with taxmatrix=NULL
data(newickstr)
z<-phyfromnewick(str=newickstr)
phyreg(phydata=z$phy,taxmatrix=NULL)
##   ... and if branch lengths were supplied, and we trust them, we can
phyreg(y~X1, "A", phydata=z$phy, heightsdata=z$hts)
##
##  Similarly with a phylogeny in phylo format. We obtain one of those by loading "exphylo".
data(exphylo)
phyconverted<-phyfromphylo(exphylo)
phyreg(phydata=phyconverted$phy, heightsdata=NULL)
## If we don't unset heightsdata, we get an error. Although phyconverted and z describe the
##   same phylogeny, there is a singleton higher node (i.e. only one daughter) in
##   phyconverted, which phyfromphylo has left in. But phyfromnewick strips them out.
##   But of course the phy and hts from the same conversion always work together:
phyreg(phydata=phyconverted$phy, heightsdata=phyconverted$hts)
## The results should all be the same (except for where we switch the main dataset) because 
##  it's the same phylogeny represented in three different ways, with the same heights.
##
## Now we obtain a plot to investigate which higher nodes are important in the influence of the 
##  test variables
m2<-phyreg(plot=1)
##    ... and we look at the data that's been plotted
m2$details$influence
##    The bottom right corner has a point that fits much worse after test is added, and the
##    table we've just printed tells us that is node 115. To find out which node that really
##    is, we just do
who(m2$fullphy)
##    and a table appears telling us which species are under each higher node. Note that
##    because phyfromphylo retains singleton higher nodes (those with just one daughter)
##    that are present in the original phylo structure, two nodes (121 and 126) have the same
##    descendant species (81:100), and one of the nodes (126) is omitted from the short
##    regression.
##
## Finally, let's do an lm() to find the true values of the parameters in the dataset myd2
lm(y~X1+X2+X3+A+B+error,data=myd2)
##   ... and all is laid bare!
##
## Enjoy!
##
## NOTA BENE: if you have actually run the examples, the parameters will have been left as 
##            they are and will carry on to your session, so you can explore the examples if
##            you like. To reset them for your own work, include "reset=TRUE" in the call
##            to phyreg.
##            
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ linear model }
\keyword{ comparative data }% __ONLY ONE__ keyword per line
