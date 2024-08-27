\name{PLOT.PCA}
\alias{PLOT.PCA}
\title{
Plot the results of a principal component analysis
}
\description{
Various maps and biplots of the results of a principal component analysis using function \code{PCA}.
}
\usage{
PLOT.PCA(obj, map="symmetric", rescale=1, dim=c(1,2), axes.inv = c(1,1),
         main="", cols=c("blue","red"), colarrows = "pink", cexs=c(0.8,0.8),
         fonts=c(2,4))
}
\arguments{
  \item{obj}{
An LRA object created using function \code{LRA}
}
 \item{map}{
Choice of scaling of rows and columns: \code{"symmetric"} (default), \code{"asymmetric"} or \code{"contribution"}
}
  \item{rescale}{
A rescaling factor applied to column coordinates (default is 1 for no rescaling)
}
 \item{dim}{
Dimensions selected for horizontal and vertical axes of the plot (default is c(1,2)) 
}
 \item{axes.inv}{
Option for reversing directions of horizontal and vertical axes (default is c(1,1) for no reversing, change one or both to -1 for reversing)
}
 \item{main}{
Title for plot
}
 \item{cols}{
Colours for row and column labels (default is c("blue","red"))
}
 \item{colarrows}{
Colour for arrows in asymmetric and contribution biplots (default is "pink")
}
 \item{cexs}{
Character expansion factors for row and column labels (default is c(0.8,0.8))
}
 \item{fonts}{
Fonts for row and column labels (default is c(2,4))
}
}
\details{
The function \code{PLOT.PCA} makes a scatterplot of the results of a logratio analysis (computed using function \code{PCA}), with various options for scaling the results and changing the direction of the axes.  By default, dimensions 1 and 2 are plotted on the horizontal and vertical axes, and it is assumed that row points refer to samples and columns to variables.

By default, the symmetric scaling is used, where both rows and columns are in principal coordinates and have the same amount of weighted variance along the two dimensions. The other options are biplots: the asymmetric option, when columns are in standard coordinates, and the contribution option, when columns are in contribution coordinates.  In cases where the row and column displays occupy widely different extents, the column coordinates can be rescaled using the \code{rescale} option. 
}
\references{
Greenacre, M. (2013), Contribution biplots, Journal of Computational and Graphical Statistics, 22, 107-122.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link[ca]{plot.ca}}
}
\examples{
# perform weighted PCA on the ALRs of the Roman cups data set
# where the first oxide silica is chosen as the denominator
data("cups")
cups.alr <- ALR(cups, denom=1)
cups.pca <- PCA(cups.alr)
PLOT.PCA(cups.pca, map="contribution", rescale=0.2, axes.inv=c(1,-1))	
}
\keyword{multivariate}
