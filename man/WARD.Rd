\name{WARD}
\alias{WARD}
\title{
Ward clustering of a compositional data matrix
}
\description{
This function clusters the rows (or the columns, if the matrix is transformed) of a compositional data matrix, using weighted Ward clustering of the logratios.
}
\usage{
WARD(LRdata, weight=TRUE, row.wt=NA)
}
\arguments{
  \item{LRdata}{
Matrix of logratios, either a vector or preferably the logratio object resulting from one of the functions ALR, CLR, PLR or LR  (usually CLRs will be used))
}
  \item{weight}{
\code{TRUE} (default) for weighted analysis (in which case weights are in the logratio object), \code{FALSE} for unweighted analysis, or a vector of user-defined column weights
}
 \item{row.wt}{
Optional set of row weights (default is equal weights when \code{row.wt=NA})
}
}
\details{
The function \code{WARD} performs a weighted WARD hierarchical clustering on the rows of an input set of logratios, usually CLR-transformed.
(This would be equivalent to performing the clustering on all pairwise logratios).
If the columns of the logratio matrix are unweighted, specify the option \code{weight=FALSE}: they will then get equal weights.
The default \code{weight=TRUE} option implies that column weights are provided, either in the input list object LRdata, as \code{LRdata$LR.wt}, or as
a vector of user-specified weights using the same \code{weight} option.
}
\value{
An object which describes the tree produced by the clustering process on the n objects. The object is a list with components:
  \item{merge  }{an n-1 by 2 matrix. Row i of \code{merge} describes the merging of clusters at step i of the clustering. If an element j in the row is negative, then observation -j was merged at this stage. If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm. Thus negative entries in \code{merge} indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.}
  \item{height }{a set of n-1 real values (non-decreasing for ultrametric trees). The clustering height: that is, the value of the criterion associated with the clustering method for the particular agglomeration.}
  \item{order  }{a vector giving the permutation of the original observations suitable for plotting, in the sense that a cluster plot using this ordering and matrix merge will not have crossings of the branches}
}
\references{
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC. \cr
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{hclust}}, \code{\link{CLR}}, \code{\link{LR.VAR}}, \code{\link{CLOSE}}
}
\examples{
# clustering steps for unweighted and weighted logratios
# for both row- and column-clustering
data(cups)
cups <- CLOSE(cups)

# unweighted logratios: clustering samples
cups.uclr  <- CLR(cups, weight=FALSE)
cups.uward <- WARD(cups.uclr, weight=FALSE)   # weight=FALSE not needed here,
                                              # as equal weights are in object
plot(cups.uward)
# add up the heights of the nodes
sum(cups.uward$height)
# [1] 0.02100676
# check against the total logratio variance
LR.VAR(cups.uclr, weight=FALSE)
# [1] 0.02100676

# unweighted logratios: clustering parts
tcups <- t(cups)
tcups.uclr  <- CLR(tcups, weight=FALSE)
tcups.uward <- WARD(tcups.uclr, weight=FALSE)  # weight=FALSE not needed here,
                                               # as equal weights are in object
plot(tcups.uward, labels=colnames(cups))
sum(tcups.uward$height)
# [1] 0.02100676
LR.VAR(tcups.uclr, weight=FALSE)
# [1] 0.02100676

# weighted logratios: clustering samples
cups.clr <- CLR(cups)
cups.ward <- WARD(cups.clr)
plot(cups.ward)
sum(cups.ward$height)
# [1] 0.002339335
LR.VAR(cups.clr)
# [1] 0.002339335

# weighted logratios: clustering parts
# weight=FALSE is needed here, since we want equal weights
# for the samples (columns of tcups)
tcups.clr  <- CLR(tcups, weight=FALSE)
tcups.ward <- WARD(tcups.clr, row.wt=colMeans(cups))
plot(tcups.ward, labels=colnames(cups))
   sum(tcups.ward$height)
# [1] 0.002339335
LR.VAR(tcups.clr, row.wt=colMeans(cups))
# [1] 0.002339335
}
\keyword{multivariate}
