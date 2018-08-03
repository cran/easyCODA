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
Matrix of logratios, either a vector or preferably the logratio object resulting from one of the functions ALR, CLR, PLR or LR
}
  \item{row.wt}{
Optional set of row weights (default: equal weights)
}
  \item{weight}{
Logical indicating if varying weights are returned(default:\code{TRUE}). If \code{FALSE}, unweighted (equal) weights are returned. Alternatively a set of positive weights can be specified.
}
  \item{vars}{
If \code{TRUE}, output individual variances as well (default \code{FALSE})
}
}
\details{
The function \code{LR.VAR} computes the sum of the logratio variances provided as input, using the weights in the logratio data object. 
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
Aitchison, J. (1986), The Statistical Analysis of Compositional Data, Chapman & Hall.  \cr
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC Press.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{LR}}, \code{\link{ALR}}, \code{\link{invALR}}, \code{\link{CLR}}, \code{\link{invCLR}}
}
\examples{
data(cups)
# These give identical total logratio variances (weighted, by default)
LR.VAR(CLR(cups))
LR.VAR(LR(cups))
# Summing over all sets of ALRs gives twice the variance
totvar <- 0
for(j in 1:ncol(cups)) totvar <- totvar + LR.VAR(ALR(cups, denom=j))
totvar/2
}
\keyword{logratio}
