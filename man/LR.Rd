\name{LR}
\alias{LR}
\title{
All pairwise logratios
}
\description{
Computation of all pairwise logratios.
}
\usage{
LR(data, ordering=1:ncol(data), weight=TRUE)
}
\arguments{
  \item{data}{
A compositional data frame or matrix 
}
  \item{ordering}{
A permutation of the columns (default: the original ordering) 
}
  \item{weight}{
Logical indicating if varying weights are returned (default:\code{TRUE}). If \code{FALSE}, unweighted (equal) weights are returned. Alternatively a set of positive weights can be specified.
}
}
\details{
The function \code{LR} computes the complete set of pairwise logratios, in the order [1,2], [1,3], [2,3], [1,4], [2,4], [3,4], etc. 
}
\value{
  \item{LR}{
  The pairwise logratios as columns of a data matrix
  }
  \item{LR.wt}{
  The weights assigned to the respective logratios
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
\code{\link{ALR}}, \code{\link{invALR}}, \code{\link{CLR}}, \code{\link{invCLR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
LR(veg)
}
\keyword{logratio}
