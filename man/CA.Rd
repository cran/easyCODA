\name{CA}
\alias{CA}
\title{
Correspondence analysis
}
\description{
Computation of correspondence analysis on a table of nonnegative data.
}
\usage{
CA(data, nd = 2, suprow = NA, supcol = NA)
}
\arguments{
  \item{data}{
A data frame or matrix of nonnegative data (no negative values)
}
 \item{nd}{
Number of dimensions for summary solution if not 2 (default)
}
 \item{suprow}{
Indices of rows that are supplementary points
}
 \item{supcol}{
Indices of columns that are supplementary points
}
}
\details{
The function \code{CA} is a simple wrapper for the \code{ca} function in the \bold{ca} package (Nenadic and Greenacre, 2007), for compatibility within the \bold{easyCODA} package.

Supplementary rows and columns can be declared (also known as passive points) -- these do not contribute to the solution but are positioned on the solution axes.

The function borrows the structure and functions of the \code{ca} package, which is required, and produces a \code{ca} object, and the same \code{print}, \code{summary} and \code{plot} methods can be used, as for a \code{ca} object.     
It additionally exports the principal coordinates of both the rows and columns, not presently found in the \code{ca} package.
}
\value{
  \item{sv        }{Singular values}
  \item{nd        }{Number of dimensions in solution results}
  \item{rownames  }{Row names}
  \item{rowmass   }{Row weights}
  \item{rowdist   }{Row logratio distances to centroid}
  \item{rowinertia}{Row inertias}
  \item{rowcoord  }{Row standard coordinates}
  \item{rowpcoord }{Row principal coordinates}
  \item{rowsup    }{Indices of row supplementary points}
  \item{colnames  }{Column names}
  \item{colmass   }{Column weights}
  \item{coldist   }{Column logratio distances to centroid}
  \item{colinertia}{Column inertias}
  \item{colcoord  }{Column standard coordinates}
  \item{colpcoord }{Column principal coordinates}
  \item{N         }{The compositional data table}
}
\references{
Nenadic, O. and Greenacre, M. (2007). Correspondence analysis in R, with two- and three-dimensional graphics: The ca package. \emph{Journal of Statistical Software}, \bold{20 (3)}, \url{http://www.jstatsoft.org/v20/i03/}
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{PLOT.CA}}, \code{\link[ca]{plot.ca}}, \code{\link[ca]{summary.ca}}, \code{\link[ca]{print.ca}}
}
\examples{
# (weighted) LRA of the RomanCups data set, showing default symmetric map
data(cups)
PLOT.CA(CA(cups))

data(author)
ca(author)
plot(CA(author))      # ca plot function
PLOT.CA(CA(author))   # easyCODA plot function
}
\keyword{multivariate}
