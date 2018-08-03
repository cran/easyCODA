\name{invALR}
\alias{invALR}
\title{
Inverse of additive logratios
}
\description{
 Given additive logratios (ALRs) with respect to a specified part, compute the inverse (i.e. original parts)
}
\usage{
invALR(ALRmatrix, part.names=paste("part",1:(ncol(ALRmatrix)+1),sep=""), denom=NA)
}
\arguments{
  \item{ALRmatrix}{
A matrix of additive logratios (ALRs) with respect to a specified part)
}
 \item{part.names}{
  Part names in the reconstructed compositional data matrix
}
  \item{denom}{
   The index of the denominator used in the computation of the ALRs (default: last part))
}
}
\details{
The function \code{invALR} computes the original parts, given the additive logratios (ALRs)
}
\value{
  \item{parts}{
  The reconstructed parts (they add up to 1)
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
\code{\link{ALR}}, \code{\link{LR}}, \code{\link{CLR}}, \code{\link{invCLR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
# compute additive logratios with respect to second part
veg.ALR <- ALR(veg, denom=2)$LR
# recover original parts (to get same order, specify the denominator used originally)
invALR(veg.ALR, denom=2)
}
\keyword{logratio}
