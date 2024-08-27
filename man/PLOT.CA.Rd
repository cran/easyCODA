\name{PLOT.CA}
\alias{PLOT.CA}
\title{
Plot the results of a correspondence analysis
}
\description{
Various maps and biplots of the results of a correspondence analysis using function \code{CA}.
}
\usage{
PLOT.CA(obj, map="symmetric", rescale=1, dim=c(1,2), axes.inv = c(1,1), main="",
        cols=c("blue","red"), colarrows = "pink", cexs=c(0.8,0.8), fonts=c(2,4))
}
\arguments{
  \item{obj}{
A CA object created using function \code{CA}
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
 \item{main}{
Title for plot
}
 \item{axes.inv}{
Option for reversing directions of horizontal and vertical axes (default is c(1,1) for no reversing, change one or both to -1 for reversing)
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
The function \code{PLOT.CA} makes a scatterplot of the results of a correspondence analysis (computed using function \code{CA}), with various options for scaling the results and changing the direction of the axes.  By default, dimensions 1 and 2 are plotted on the horizontal and vertical axes, and it is assumed that row points refer to samples and columns to variables.

By default, the symmetric scaling is used, where both rows and columns are in principal coordinates and have the same amount of weighted variance (i.e. inertia) along the two dimensions. The other options are biplots: the asymmetric option, when columns are in standard coordinates, and the contribution option, when columns are in contribution coordinates.  In cases where the row and column displays occupy widely different extents, the column coordinates can be rescaled using the \code{rescale} option. 
}
\references{
Greenacre, M. (2013), Contribution biplots, Journal of Computational and Graphical Statistics, 22, 107-122.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{CA}}, \code{\link[ca]{plot.ca}}
}
\examples{
data("cups")
cups.ca <- CA(cups)
PLOT.CA(cups.ca, map="contribution", rescale=0.2)	
# Compare the above plot with that of a weighted LRA -- practically the same
cups.lra <- LRA(cups)
PLOT.LRA(cups.lra, map="contribution", rescale=0.2)
}
\keyword{multivariate}
