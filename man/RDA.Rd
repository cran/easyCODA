\name{RDA}
\alias{RDA}
\title{
Redundancy analysis
}
\description{
Computation of weighted or unweighted redundancy analysis of a samples-by-parts compositional data table, given a set of covariates.
}
\usage{
RDA(data, cov=NA, nd = NA, weight = TRUE, suprow = NA, row.wt = NA) 
}
\arguments{
  \item{data}{
A data frame or matrix of interval-scale data, e.g. logratios (which are preferably in a list object with weights)
}
 \item{cov}{
List of covariates for constraining solution
}
 \item{nd}{
Number of dimensions for summary output, by default the number of constraining dimensions
}
  \item{weight}{
TRUE (default) when weights are in data list object, FALSE for unweighted analysis, or a vector of user-defined part weights
}
 \item{suprow}{
Indices of rows that are supplementary (passive) points  (NOTE: this option is not implemented in this version)
}
 \item{row.wt}{
Optional user-defined set of positive weights for the rows (samples) (default: equal weights)
}
}
\details{
The function \code{RDA} computes a redundancy analysis of a matrix of interval-scaled data, constrained by a matrix of covariates, using the singular value decomposition. 
By default weights are assumed in the data list object. 
For the unweighted logratio analysis, specify the option \code{weight=FALSE}. 
If \code{weight = TRUE} (the default) it is assumed that the weights are included in the data object, which comes from one of the logratio functions.
User-specified part weights can be provided using the same \code{weight} option. 

Usually row weights are not specified, they are equal unless intentional weighting of the samples is desired.
Supplementary rows can be declared (also known as passive points) -- these do not contribute to the solution but are positioned on the solution axes.
This option will be available in the next release of the package.
}
\value{
  \item{sv        }{Singular values}
  \item{nd        }{Number of dimensions in the solution output}
  \item{rownames  }{Row names}
  \item{rowmass   }{Row weights}
  \item{rowdist   }{Row distances to centroid}
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
  \item{covcoord  }{Regression coordinates of constraining variables}
  \item{covnames  }{Names of constraining variables}
  \item{N         }{The data table (usually logratios in this package)}
  \item{cov       }{The table of covariates}
}
\references{
Van den Wollenbergh, A. (1977), Redundancy analysis. An alternative to canonical correlation analysis, Psychometrika 42, 207-219.\cr
Greenacre, M. (2013), Contribution biplots, Journal of Computational and Graphical Statistics 22, 107-122.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{PLOT.RDA}}, \code{\link{CLR}}, \code{\link{LR}}, \code{\link{DUMMY}}
}
\examples{
# Data frame fish has sex, habitat and mass in first columns, 
# then morphometric data in remaining columns
data(fish)
sex     <- fish[,1]
habitat <- fish[,2]
mass    <- fish[,3]
fishm   <- as.matrix(fish[,4:29])
# Convert to compositional data matrix
fishm   <- fishm / apply(fishm, 1, sum)
# Compute logarithm of mass and interaction of sex (F/M) and habitat (L/P) categories
logmass <- log(mass)
sexhab  <- 2*(sex-1)+habitat
sexhab.names <- c("FL","FP","ML","MP")
rownames(fishm) <- sexhab.names[sexhab]
# Create dummy variables for sexhab and create matrix of covariates
sexhab.Z <- DUMMY(sexhab, catnames=sexhab.names)
vars     <- cbind(logmass, sexhab.Z)
# Perform RDA on centred logratios
require(ca)
fish.RDA <- RDA(CLR(fishm), cov=vars)
# Plot results 
# (for more options see Appendix of Compositional Data Analysis in Practice)
PLOT.RDA(fish.RDA, map="contribution", rescale=0.05, indcat=2:5, 
         colrows=rainbow(4, start=0.1, end=0.8)[sexhab], cexs=c(0.8,0.8,1))
}
\keyword{multivariate}
