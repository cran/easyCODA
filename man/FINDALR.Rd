\name{FINDALR}
\alias{FINDALR}
\title{
Find the best ALR transformation
}
\description{
Searching over every possible reference part for choosing an optimal ALR transformation.
}
\usage{
FINDALR(data, weight=FALSE)
}
\arguments{
  \item{data}{
Compositional data matrix, with the parts as columns
}
  \item{weight}{
FALSE (default) for equally weighted parts, TRUE when weights are in data list object, or a vector of user-defined part weights
}
}
\details{
The function \code{FINDALR} considers every possible set of additive logratio (ALR) transformations, by trying each of the references. For each set the closeness to isometry is measured by the Procrustes correlation. In addition, the variance of the log-transformed reference is also computed. The reference with highest Procrustes correlation and the reference with the lowest variance of its log-transform are identified. The number of ALRs computed is equal to 1 less than the number of rows or columns, whichever is the smallest.
}
\value{
An object which describes the tree produced by the clustering process on the n objects. The object is a list with components:
  \item{totvar}{Total logratio variance}
  \item{procrust.cor }{The Procrustes correlations of the ALRs using each reference}
  \item{procrust.max}{The value of the highest Procrustes correlation}
  \item{procrust.ref}{The reference corresponding to the highest correlation}
  \item{var.log}{Variances of the log-transformed references}
  \item{var.min}{The value of the lowest variance}
  \item{var.ref}{The reference corresponding to the lowest variance}
}
\references{
Greenacre, M., Martinez-Alvaro, M. and Blasco, A. (2021), Compositional data analysis of microbiome and any-omics datasets: a validation of the additive logratio transformation, Frontiers in Microbiology 12: 2625 \cr
Gower, J. and Dijksterhuis, G.B. (2004), Procrustes Problems. Oxford University Press   \cr
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{ALR}}
}
\examples{
# For the fish morphometric data, first close (normalize) 
# then loop over the 26 possible references
data(fish)
FINDALR(CLOSE(fish[,4:29]))
# Note that for the default option weight=FALSE closing the data is not necessary
}
\keyword{multivariate}
