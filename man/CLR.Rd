\name{CLR}
\alias{CLR}
\title{
Centred logratios
}
\description{
Computation of centred logratios (CLRs).
}
\usage{
CLR(data, weight=TRUE)
}
\arguments{
  \item{data}{
A compositional data frame or matrix 
}
  \item{weight}{
Logical indicating if varying weights are returned(default:\code{TRUE}). If \code{FALSE}, unweighted (equal) weights are returned. Alternatively a set of positive weights can be specified.
}
}
\details{
The function \code{CLR} computes the set of centred logratios (CLRs). 
}
\value{
  \item{LR}{
  The centred logratios (CLRs)
  }
  \item{LR.wt}{
  The weights assigned to the CLRs
  }
}
\references{
Aitchison, J. (1986), The Statistical Analysis of Compositional Data, Chapman & Hall. \cr
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC Press.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{invCLR}}, \code{\link{ALR}}, \code{\link{invALR}}, \code{\link{LR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
CLR(veg)
}
\keyword{logratio}
