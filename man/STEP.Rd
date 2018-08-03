\name{STEP}
\alias{STEP}
\title{
Stepwise selection of logratios
}
\description{
Stepwise selection of pairwise logratios that explain maximum variance in a target matrix.
}
\usage{
STEP(data, datatarget=data, previous=NA, previous.wt=NA, weight=TRUE, 
     random=FALSE, nsteps=min(ncol(data), ncol(datatarget))-1, top=1)
}
\arguments{
  \item{data}{
A data frame or matrix of compositional data on which pairwise logratios are computed
}
 \item{datatarget}{
A matrix of interval-scale data, with as many rows as \code{data}, which serves as the target matrix whose variance is to be explained (by default it is the same matrix as data, in which case total logratio variance is to be explained)
}
 \item{previous}{
A vector or matrix of variables to be forced in before logratios are sought
}
 \item{previous.wt}{
Possible weights of the variable(s) forced in before logratios are sought (if not specified, weights of 1 are assumed)
}
  \item{weight}{
\code{TRUE} (default) when weights are in data list object, \code{FALSE} for unweighted analysis, or a vector of user-defined part weights
}
 \item{random}{
\code{TRUE} if a random selection is made of tied logratios; \code{FALSE} (default) if logratio that maximizes Procrustes correlation is chosen  
}
 \item{nsteps}{
Number of steps to take (by default, one less than the number of columns of data and of datatarget, whichever is smaller)
}
 \item{top}{
Number of top variance-explaining logratios returned after last step (by default, 1, i.e. the best) 
}
}
\details{
The function \code{STEP} sequentially computes the logratios in a data matrix (usually compositional) that best explain the variance in a second matrix, called the target matrix.  By default, the target matrix is the same matrix, in which case the logratios that best explain the logratio variance in the same matrix are computed.
In this case, weights for the data matrix are assumed by default, proportional to part means of the compositional data matrix.
For the unweighted logratio variance, specify the option \code{weight=FALSE}. 
User-specified weights on the columns of the data matrix (usually compositional parts) can be provided using the same \code{weight} option. 

If the target matrix is a different matrix, it is assumed that its regular total variance is to be explained, not its logratio variance.

If \code{nsteps > 1} and \code{top=1} the results are in the form of an optimal set of logratios that sequentially add maximum explained variance at each step.  If \code{top>1} then at the last step the ordered list of top variance-explaining logratios is returned, which allows users to make an alternative choice of the logratio based on substantive knowledge. Hence, if \code{nsteps=1} and \code{top=10}, for example, the procedure will move only one step, but list the top 10 logratios for that step.
}
\value{
  \item{rationames }{Names of maximizing ratios in stepwise process}
  \item{ratios     }{Indices of ratios}
  \item{logratios  }{Matrix of logratios}
  \item{R2max      }{Sequence of maximum cumulative explained variances}
  \item{procr      }{Corresponding sequence of Procrustes correlations}
  \item{rationames.top }{Names of "top" ratios at last step}
  \item{ratios.top     }{Indices of "top" ratios}
  \item{logratios.top  }{Matrix of "top" logratios}
  \item{R2.top         }{Sequence of "top" cumulative explained variances (in descending order)}
  \item{procr.top      }{Corresponding sequence of "top" Procrustes correlations}
}
\references{
Van den Wollenbergh, A. (1977), Redundancy analysis. An alternative to canonical correlation analysis, Psychometrika 42, 207-219.\cr
Greenacre, M. (2018), Variable selection in compositional data analysis using pairwise logratios, Mathematical Geosciences, DOI: 10.1007/s11004-018-9754-x. \cr
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{PLOT.RDA}}, \code{\link{CLR}}, \code{\link{LR}}, \code{\link{ALR}}
}
\examples{
# Stepwise selection of ratios for RomanCups data set
data(cups)
# Set seed to obtain same results as in Appendix C of Greenacre (2018)
set.seed(2872)
STEP(cups, random=TRUE)
# Select best ratio, but output "top 5"
STEP(cups, nsteps=1, top=5)
}
\keyword{multivariate}
