\name{DOT}
\alias{DOT}
\title{
Dot plot
}
\description{
Simple dot plot of original data
}
\usage{
DOT(data, cols=NA, names=NA, groups=NA, pch=19, horizon=FALSE, jitter=1, 
     xscale=NA, xscalefac=1, yaxis=TRUE, shownames=TRUE, main="", ylab="", 
     xlim=c(NA,NA), ylim=c(NA, NA), cex=1)
}
\arguments{
  \item{data}{
Matrix or data frame with data in columns (even if single sample, should be a matrix)
}
 \item{cols}{
Colours of points for each sample, default rainbow
}
  \item{names}{
Labels for variables, by default the column names of data, or group names
}
 \item{groups}{
Group codes to split the data vector into separate plots
}
 \item{pch}{
Point character
}
 \item{horizon}{
\code{TRUE} if horizontal gray dashed lines required at "nice" y-values (default FALSE, not implemented yet)
}
 \item{jitter}{
1 by default, increase or decrease slightly for more jitter
}
 \item{xscale}{
User-supplied positions of points on horizontal axis
}
 \item{xscalefac}{
1 by default, rescale the positions on horizontal axis
}
 \item{yaxis}{
TRUE by default, FALSE to suppress and optionally add afterwards
}
  \item{shownames}{
\code{TRUE} by default; \code{FALSE} to not show group names and add them externally
}
  \item{main}{
Heading
}
  \item{ylab}{
Vertical axis label
}
  \item{xlim}{
Horizontal axis limits 
}
  \item{ylim}{
Vertical axis limits 
}
 \item{cex}{
Character size adjustment for labels
}
}
\details{
The function \code{DOT} makes a dot plot for specified groups of points, which can be in columns of a matrix or data frame. 
}
\references{
Greenacre, M. (2016), Data reporting and visualization in ecology, Polar Biology, 39:2189-2205.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{BAR}}
}
\examples{
data(veg)
DOT(veg)
}
\keyword{visualization}
