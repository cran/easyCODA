\name{ACLUST}
\alias{ACLUST}
\title{
Amalgamation clustering of the parts of a compositional data matrix
}
\description{
This function clusters the parts of a compositional data matrix, using amalgamation of the parts at each step.
}
\usage{
ACLUST(data, weight = TRUE, close = TRUE)
}
\arguments{
  \item{data}{
Compositional data matrix, with the parts as columns
}
  \item{weight}{
\code{TRUE} (default) for weighting using part averages of closed compositions, \code{FALSE} for unweighted analysis, or a vector of user-defined column weights
}
  \item{close}{
\code{TRUE} (default) will close the rows of \code{data} prior to clustering, \code{FALSE} leaves \code{data} as it is
}
}
\details{
The function \code{ACLUST} performs amalgamation hierarchical clustering on the parts (columns) of a given compositional data matrix, as proposed by Greenacre (2019).
At each step of the clustering two clusters are amalgamated that give the least loss of explained logratio variance.
}
\value{
An object which describes the tree produced by the clustering process on the n objects. The object is a list with components:
  \item{merge  }{an n-1 by 2 matrix. Row i of \code{merge} describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in \code{merge} indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.}
  \item{height }{a set of n-1 real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.}
  \item{order  }{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches}
  \item{labels }{a vector of column labels, the column names of \code{data}}
}
\references{
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC. \cr
Greenacre, M. (2019), Amalgamations are valid in compositional data analysis,
                      can be used in agglomerative clustering, and their
                      logratios have an inverse transformation.
                      Applied Computing and Geosciences, open access.\cr
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{hclust}}, \code{\link{WARD}},\code{\link{CLR}}, \code{\link{LR.VAR}},
\code{\link{CLOSE}}
}
\examples{
data(cups)

# amalgamation clustering    (weighted parts)
cups.aclust <- ACLUST(cups)
plot(cups.aclust)

# reproducing Figure 2(b) of Greenacre (2019) (unweighted parts))
# dataset Aar is in the compositions package
# aar is a subset of Aar
# code given here within the '\dontrun' environment since external package 'compositions' required
\dontrun{
  library(compositions)
  data(Aar)
  aar <- Aar[,c(3:12)]
  aar.aclust <- ACLUST(aar, weight=FALSE)
# the maximum height is the total variance
# convert to percents of variance NOT explained
  aar.aclust$height <- 100 * aar.aclust$height / max(aar.aclust$height)
  plot(aar.aclust, main="Parts of Unexplained Variance", ylab="Variance (percent)")
}
}
\keyword{multivariate}
