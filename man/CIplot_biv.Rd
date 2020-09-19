\name{CIplot_biv}
\alias{CIplot_biv}
\title{
Bivariate confidence and data ellipses
}
\description{
Draws confidence and data ellipses in bivariate scatterplots
}
\usage{
CIplot_biv(x, y, group, wt=rep(1/length(x),length(x)), 
           varnames=c("x","y"), groupnames=sort(unique(group)), 
           groupcols=rainbow(length(unique(group))), 
           shownames=TRUE, xlim=c(NA,NA), ylim=c(NA,NA), 
           lty=1, lwd=1, add=FALSE, alpha=0.95, ellipse=0, 
           shade=FALSE, frac=0.01, cex=1) 
}
\arguments{
  \item{x}{
x-variable (horizontal) of scatterplot
}
 \item{y}{
y-variable (vertical) of scatterplot
}
  \item{group}{
Grouping variable
}
 \item{wt}{
Set of weights on the cases (operates when ellipse=1)
}
 \item{varnames}{
Vector of two labels for the axes (default is x and y)
}
 \item{groupnames}{
Vector of labels for the groups (default is 1, 2, etc...)
}
 \item{groupcols}{
Vector of colours for the groups
}
 \item{shownames}{
Whether to show group names at group centroids or not (default is \code{TRUE})
}
 \item{xlim}{
Possible new x-limits for plot
}
 \item{ylim}{
Possible new y-limits for plot
}
  \item{lty}{
Line type for the ellipses (default is 1)
}
 \item{lwd}{
Line width for the ellipses (default is 1)
}
  \item{add}{
=\code{TRUE} if ellipses/intervals are added to existing plot (default is \code{FALSE})
}
  \item{alpha}{
Confidence level of ellipses (default is 0.95)
}
  \item{ellipse}{
Type of ellipse (see Details below; default is 0 for normal-based ellipses)
}
  \item{shade}{
=\code{TRUE} for ellipse shading (default=\code{FALSE})
}
  \item{frac}{
Proportional part defining the width of the bars at the edges of confidence intervals (for \code{ellipse=3} and \code{4})
}
 \item{cex}{
Character expansion factor for group names
}
}
\details{
The function \code{CIplot_biv} makes various types of confidence and data ellipses, according to option \code{ellipse}. Set \code{ellipse<0} for regular data-covering ellipses. Set \code{ellipse=0} (default) for normal-theory confidence ellipses. Set\code{ellipse=1} for bootstrap confidence ellipses. The option \code{ellipse=2} for the delta method is not implemented yet.  Set \code{ellipse=3} for normal-theory confidence error bars lined up with axes. Set \code{ellipse=4} for bootstrap confidence error bars along axes.  The package \code{ellipse} is required.

}
\references{
Greenacre, M. (2016), Data reporting and visualization in ecology, Polar Biology, 39:2189-2205.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{DOT}}
}
\examples{
# Generate some bivariate normal data in three groups with different means
# Means (1,0), (0,1) and (0,0)
means  <- matrix(c(1,0,0,1,0,0), ncol=3)
data   <- matrix(nrow=300, ncol=2)
groups <- sample(rep(c(1,2,3), 100))
for(i in 1:300) data[i,] <- rnorm(c(1,1), mean=means[,groups[i]])
# Plot confidence ellipses with shading
CIplot_biv(data[,1], data[,2], group=groups, shade=TRUE)
}
\keyword{visualization}
