\name{ILR}
\alias{ILR}
\title{
Isometric logratio
}
\description{
Computation of a single isometric logratio (ILR)
}
\usage{
ILR(data, numer=NA, denom=NA, weight=TRUE)
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
The function \code{ILR} computes a single isometric logratio based on the specified numerator and denominator parts that define the two geometric means in the ratio. 
}
\value{
  \item{LR}{
  The isometric logratio (ILR)
  }
  \item{LR.wt}{
  The weight assigned to the ILR
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
\code{\link{SLR}}, \code{\link{ALR}}, \code{\link{PLR}}, \code{\link{LR}}
}
\examples{
data(veg)
ILR(veg, numer=1, denom=2:3)
}
\keyword{logratio}
