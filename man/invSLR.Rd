\name{invSLR}
\alias{invSLR}
\title{
Inverse of full set of amalgamation balances
}
\description{
Given a full set of amalgamation (or summation) balances (SLRs), compute the inverse (i.e. recover the original parts)
}
\usage{
invSLR(SLRmatrix, part.names=NA, ratio.names=colnames(SLRmatrix))
}
\arguments{
  \item{SLRmatrix}{
  A matrix of amalgamation logratios, one less column than the number of parts
}
 \item{part.names}{
  Part names in the reconstructed compositional data matrix
}
 \item{ratio.names}{
  Definition of the amalgamation logratios
}
}
\details{
The function \code{invSLR} computes the original parts, given the amalgamation logratios (CLRs).
The amalgamation logratios are specified in \code{ratio.names} in the format \code{num/den} where \code{num} and \code{den} are the numerator and denominator amalgamations respectively. 
An amalgamation is specified as \code{"p1&p2&..."}, where \code{p1}, \code{p2}, etc. are the parts summed in the amalgamation.
For example, an SLR of the ratio MnO/(CaO+P2O5) would be names as \code{"MnO/CaO&P2O5"}. 
}
\value{
  \item{parts}{
  The reconstructed parts (they add up to 1)
  }
}
\references{
Aitchison, J. (1986), The Statistical Analysis of Compositional Data, Chapman & Hall.   \cr
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC Press.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{CLR}}, \code{\link{ALR}}, \code{\link{invALR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
# compute centred logratios
veg.CLR <- CLR(veg)$LR
# invert back to original parts (parts closed to sum to 1)
invALR(veg.CLR)
}
\keyword{logratio}
