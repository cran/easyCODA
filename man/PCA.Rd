\name{PCA}
\alias{PCA}
\title{
Principal component analysis
}
\description{
Computation of weighted or unweighted principal component analysis of a matrix of interval-scale data (e.g. a matrix of logratios).
}
\usage{
PCA(data, nd = 2, weight = TRUE, row.wt = NA, suprow = NA, supcol = NA)
}
\arguments{
  \item{data}{
A data frame or matrix of interval-scale data, or logratio object from functions \code{ALR}, \code{CLR} or \code{LR}
}
 \item{nd}{
Number of dimensions for summary solution if not 2
}
  \item{weight}{
TRUE (default) for column weighting, FALSE for unweighted analysis, or a vector of user-defined column weights
}
 \item{row.wt}{
Optional user-defined set of positive weights for the rows (samples) (default: equal weights)
}
 \item{suprow}{
Indices of rows that are supplementary points (NOTE: this option is not implemented in this version)
}
 \item{supcol}{
Indices of columns that are supplementary points (NOTE: this option is not implemented in this version)
}
}
\details{
The function \code{PCA} computes an unstandardized principal component analysis, based on the singular value decomposition, of a matrix of interval-scale data, usually a matrix of logratios in the context of this package (but it can be used for general data as well). 
For general usage the unweighted option \code{weight = FALSE} might be preferred, but the default is weighted in the present context of compositional data.

User-specified weights can be supplied, for the rows and/or the columns.
Usually row weights are not specified, and are equal unless intentional weighting of the samples is desired.
User-specified part weights can be provided using the \code{weight} option. 

Supplementary rows and columns can be declared (also known as passive points) -- these do not contribute to the solution but are positioned on the solution axes. Notice that this optyion is not implemented in the present version, but will appear in the next one.

The function borrows the structure and functions of the \code{ca} package, which is required, and produces a \code{ca} object, and the same \code{print}, \code{summary} and \code{plot} methods can be used, as for a \code{ca} object.}
\value{
  \item{sv        }{Singular values}
  \item{nd        }{Dimenson of the solution}
  \item{rownames  }{Row names}
  \item{rowmass   }{Row weights}
  \item{rowdist   }{Row logratio distances to centroid}
  \item{rowinertia}{Row variances}
  \item{rowcoord  }{Row standard coordinates}
  \item{rowpcoord }{Row principal coordinates}
  \item{rowsup    }{Indices of row supplementary points}
  \item{colnames  }{Column names}
  \item{colmass   }{Column weights}
  \item{coldist   }{Column logratio distances to centroid}
  \item{colinertia}{Column variances}
  \item{colcoord  }{Column standard coordinates}
  \item{colpcoord }{Column principal coordinates}
  \item{N         }{The data table}
}
\references{
Aitchison, J. and Greenacre, M. (2002), Biplots of compositional data, Applied Statistics 51, 375-392.\cr
Greenacre, M. (2010), Biplots in Practice, BBVA Foundation, Bilbao. Free download from www.multivariatestatistics.org
}
\author{
Michael Greenacre
}
\seealso{
, \code{PLOT.PCA}, \code{\link[ca]{plot.ca}}, \code{\link[ca]{summary.ca}}, \code{\link[ca]{print.ca}}
}
\examples{
# compute logratios of Vegetables data set
data(veg)
veg.LR <- LR(veg)
# unweighted PCA biplot of the results
veg.PCA <- PCA(veg.LR$LR, weight=FALSE)
PLOT.PCA(veg.PCA, map="asymmetric")
}
\keyword{multivariate}
