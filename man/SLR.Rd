\name{SLR}
\alias{SLR}
\title{
Amalgamation (summed) logratio
}
\description{
Computation of a single amalgamation (summed) logratio
}
\usage{
SLR(data, numer=NA, denom=NA, weight=TRUE)
}
\arguments{
  \item{data}{
A compositional data frame or matrix 
}
  \item{numer}{
Vector of parts in the numerator
}
  \item{denom}{
Vector of parts in the denominator
}
  \item{weight}{
Logical indicating if a varying weight is returned (default:\code{TRUE}). If \code{FALSE}, a weight based on equally-weighted parts is returned. Alternatively a positive weight can be specified.
}
}
\details{
The function \code{SLR} computes a single amalgamation logratio based on the specified numerator and denominator parts that define the two summations in the ratio. 
}
\value{
  \item{LR}{
  The amalgamation (summed)) logratio (SLR)
  }
  \item{LR.wt}{
  The weight assigned to the SLR
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
\code{\link{ILR}}, \code{\link{ALR}}, \code{\link{CLR}}, \code{\link{PLR}}, \code{\link{LR}}
}
\examples{
data(veg)
SLR(veg, numer=1, denom=2:3)
}
\keyword{logratio}
