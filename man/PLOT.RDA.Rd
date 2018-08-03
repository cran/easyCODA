\name{PLOT.RDA}
\alias{PLOT.RDA}
\title{
Plot the results of a redundancy analysis
}
\description{
Various maps and biplots/triplots of the results of a redundancy analysis using function \code{RDA}.
}
\usage{
PLOT.RDA(obj, map="symmetric", indcat=NA, rescale=1, dim=c(1,2), main="", 
         axes.inv=c(1,1), rowstyle=1, cols=c("blue","red","forestgreen"), 
         colarrows=c("pink","lightgreen"), colrows=NA, pchrows=NA, colcat=NA,
         cexs=c(0.8,0.8,0.8), fonts=c(2,4,4))
}
\arguments{
  \item{obj}{
An RDA object created using function \code{RDA}
}
 \item{map}{
Choice of scaling of rows and columns: \code{"symmetric"} (default), \code{"asymmetric"} or \code{"contribution"}
}
  \item{indcat}{
A vector indicating which of the covariates are dummy (or fuzzy) variables
}
  \item{rescale}{
A rescaling factor applied to column coordinates(default is 1 for no rescaling). If rescale is a vector with two values, the first applies to the column coordinates and the second to the covariate coordinates.
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
 \item{rowstyle}{
Scaling option for row coordinates, either 1 (SVD coordinates, default) or 2 (as supplementary points)
}
 \item{cols}{
Colours for row and column  and covariate labels (default is c("blue","red","forestgreen"))
}
 \item{colarrows}{
Colour for arrows in asymmetric or contribution biplots, for columns and covariates (default is c("pink","lightgreen"))
}
 \item{colrows}{
Optional vector of colours for rows
}
 \item{pchrows}{
Optional vector of point symbols for rows
}
 \item{colcat}{
Optional vector of colours for covariate categories (dummy variables)
}
 \item{cexs}{
Vector of character expansion factors for row and column and covariate labels (default is c(0.8,0.8,0.8))
}
 \item{fonts}{
Vector of font styles for row and column and covariate labels (default is c(2,4,4))
}
}
\details{
The function \code{PLOT.RDA} makes a scatterplot of the results of a redundancy analysis (computed using function \code{RDA}), with various options for scaling the results and changing the direction of the axes.  By default, dimensions 1 and 2 are plotted on the horizontal and vertical axes, and it is assumed that row points refer to samples and columns to compositional parts.   Covariates are plotted according to their regression coefficients with the RDA dimensions, and if they contain dummy (or fuzzy) variables these are indicated by the option indcat, and hence plotted as centroids not arrows.

By default, the symmetric scaling is used, where both rows and columns are in principal coordinates and have the same amount of weighted variance along the two dimensions. The other options are the asymmetric option, when columns are in standard coordinates, and the contribution option, when columns are in contribution coordinates.  In cases where the row and column displays as well as the covariate positions occupy widely different extents, the column and covariate coordinates can be rescaled using the \code{rescale} option. 
}

\author{
Michael Greenacre
}
\seealso{
\code{\link{RDA}}
}
\examples{
# see the use of PLOT.RDA in the example of the RDA function
}
\keyword{multivariate}
