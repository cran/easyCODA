\name{ALR}
\alias{ALR}
\title{
Additive logratios
}
\description{
Computation of additive logratios (ALRs) with respect to a specified part.
}
\usage{
ALR(data, denom=ncol(data), weight=TRUE, stats=FALSE)
}
\arguments{
  \item{data}{
A compositional data frame or matrix 
}
 \item{denom}{
Number of part used in the denominator
}
  \item{weight}{
Logical indicating if varying weights are returned(default:\code{TRUE}). If \code{FALSE}, unweighted (equal) weights are returned. Alternatively a set of positive weights can be specified.
}
 \item{stats}{
Logical indicating if means, variances and total variance of the ALRs are returned (default:\code{FALSE})
}
}
\details{
The function \code{ALR} computes a set of additive logratios (ALRs) with respect to a specified part (by default, the last part). 
}
\value{
  \item{LR}{
  The additive logratios (ALRs)
  }
  \item{LR.wt}{
  The weights assigned to the ALRs
  }
  \item{denom}{
  The index of the denominator used in the computation of the ALRs 
  }
  \item{part.names}{
  The part names in the data, i.e. column names 
  }
  \item{part.wt}{
  The part weights 
  }
  \item{means}{
  The means of the ALRs (only returned if \code{stats = TRUE})
  }
  \item{vars}{
  The variances of the ALRs (only returned if \code{stats = TRUE})
  }
  \item{totvar}{
  The total variance of the ALRs (only returned if \code{stats = TRUE})
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
\code{\link{invALR}}, \code{\link{LR}}, \code{\link{CLR}}, \code{\link{invCLR}}, \code{\link{LR.VAR}}
}
\examples{
data(veg)
ALR(veg, denom=2)
}
\keyword{logratio}
