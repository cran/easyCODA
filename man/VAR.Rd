\name{VAR}
\alias{VAR}
\title{
Variance of a vector of observations, dividing by n rather than n-1
}
\description{
This function computes the usual variance but divides by n, not by n-1.
}
\usage{
VAR(x)
}
\arguments{
  \item{x}{
Vector of values for which variance is computed
}
}
\details{
To think of each of n observations weighted by 1/n this function \code{VAR} computes squared deviations from the mean and averages them.
Thus, the sum of squared deviations is divided by n rather than by n-1, as for the unbiased estimate of the variance.
}
\value{
The value of the variance.
}
\references{
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC. \cr
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{LR.VAR}}, \code{\link{CLOSE}}
}
\examples{
data(cups)
cups <- CLOSE(cups)

# variances using base R function var
apply(cups, 2, var)

# variances using easyCODA function VAR
apply(cups, 2, VAR)
}

\keyword{variance}
