\name{PLR}
\alias{PLR}
\title{
Pivot logratios
}
\description{
Computation of the set of pivot logratios(PLRs) based on the specified ordering of parts
}
\usage{
PLR(data, ordering=1:ncol(data), weight=TRUE)
}
\arguments{
  \item{data}{
A compositional data frame or matrix 
}
  \item{ordering}{
The ordering of the parts to be used in the PLRs (by default, the original ordering of the columns)
}
  \item{weight}{
Logical indicating if varying weights are returned (default:\code{TRUE}). If \code{FALSE}, weights based on equally-weighted parts are returned. Alternatively a set of positive weights can be specified.
}
}
\details{
The function \code{PLR} computes the set of pivot logratios according to the ordering of the parts. 
}
\value{
  \item{LR}{
  The pivot logratios (PLRs)
  }
  \item{LR.wt}{
  The weights assigned to the PLRs
  }
}
\references{
Hron K., Filzmoser P., de Caritat P., Fiserova E., Gardlo A. (2017). Weighted pivot coordinates for copositional data and their application to geochemical mapping. Mathematical Geosciences 49, 777-796.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{ILR}}, \code{\link{SLR}}, \code{\link{CLR}}, \code{\link{ALR}}, \code{\link{LR}}
}                   
\examples{
data(veg)
PLR(veg, ordering=c(1,3,2))
}
\keyword{logratio}
