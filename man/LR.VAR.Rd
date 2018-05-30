\name{LR.VAR}
\alias{LR.VAR}
\title{
Total logratio variance
}
\description{
Computation of total (weighted)logratio variance.
}
\usage{
LR.VAR(LRdata, row.wt = NA, weight=TRUE, vars=FALSE)
}
\arguments{
  \item{LRdata}{
Matrix of logratios, either a vector or preferably the output from one of the functions ALR, CLR, ILR, SLR, PLR or LR
}
  \item{row.wt}{
Optional set of row weights (default: equal weights)
}
  \item{weight}{
Logical indicating if varying weights are returned(default:\code{TRUE}). If \code{FALSE} unweighted (equal) weights are returned. Alternatively a set of positive weights can be specified.
}
  \item{vars}{
Output individual variances as well (default: no output of individual variances)
}
}
\details{
The function \code{LR} computes the comoplete set of pairwise logratios, in order ([1,2], [1,3], [2,3], [1,4], [2,4], [3,4], etc). 
}
\value{
  \item{LRtotvar}{
  The total logratio variance
  }
  \item{LRvars}{
  (optional, if vars=TRUE, the individual logratio variances composing the total)
  }
}
\references{
Aitchison, J. (1986), The Statistical Analysis of Compositional Data, Chapman & Hall. 
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC Press.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{LR}}, \code{\link{ALR}}, \code{\link{invALR}}, \code{\link{CLR}}, \code{\link{invCLR}}
}
\examples{
data(veg)
LR(veg)
}
\keyword{logratio}
