\name{BAR}
\alias{BAR}
\title{
Compositional bar plot
}
\description{
Horizontal bar plot of compositional data
}
\usage{
BAR(data, cols=rainbow(ncol(data)), col.names=colnames(data),
    row.names=rownames(data), order.column=NA, eps=0.5, main="",
    ylab="", ylim=c(0,nrow(data)), xlim=c(0,100), cex=1)
}
\arguments{
  \item{data}{
Compositional data matrix or data frame with compositions in rows, parts in columns
}
 \item{cols}{
Colours of points for each sample, default rainbow
}
 \item{col.names}{
Part names, if modified
}
  \item{row.names}{
Sample names, if modified
}
 \item{order.column}{
By default parts are taken in order of columns, but can be re-ordered using this option
}
 \item{eps}{
Small space between bars, can be modified
}
  \item{main}{
Heading
}
  \item{ylab}{
Vertical axis label
}
  \item{ylim}{
Vertical axis limits (default is the number of rows in data)
}
  \item{xlim}{
Horizontal axis limits (default c(0,100))
}
  \item{cex}{
Character size scaling factor for labels
}
}
\details{
The function \code{BAR} makes a BAR plot for specified groups of points, which can be in columns of a matrix or data frame. 
}
\references{
Greenacre, M. (2016), Data reporting and visualization in ecology, Polar Biology: 39, 2189-2205.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{DOT}}
}
\examples{
data(veg)
BAR(veg)
}
\keyword{visualization}
