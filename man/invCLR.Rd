\name{invCLR}
\alias{invCLR}
\title{
Inverse of centred logratios
}
\description{
 Given centred logratios (CLRs), compute the inverse (i.e. recover the original parts)
}
\usage{
invCLR(CLRmatrix, part.names=colnames(CLRmatrix))
}
\arguments{
  \item{CLRmatrix}{
A matrix of centred logratios
}
 \item{part.names}{
  Part names in the reconstructed compositional data matrix
}
}
\details{
The function \code{invCLR} computes the original parts, given the centred logratios (CLRs)
}
\value{
  \item{parts}{
  The reconstructed parts (they add up to 1)
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
\code{\link{CLR}}, \code{\link{invCLR}}, \code{\link{ALR}}, \code{\link{invALR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
# compute centred logratios
veg.CLR <- CLR(veg)$LR
# invert back to original parts (parts closed to sum to 1)
invALR(veg.CLR)
}
\keyword{logratio}
