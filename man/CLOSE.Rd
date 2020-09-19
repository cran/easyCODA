\name{CLOSE}
\alias{CLOSE}
\title{
Closure of rows of compositional data matrix
}
\description{
This function closes (or normalizes) the rows of a compositional data matrix, resulting in rows summing to 1.
}
\usage{
CLOSE(x)
}
\arguments{
  \item{x}{
Compositional data matrix.
}
}
\details{
Compositional data carry relative information. It is sometimes required to close the data so that each row of observations sums to 1.
The function \code{CLOSE} performs the closure.
}
\value{
The closed compositional data matrix.
}
\references{
Greenacre, M. (2018), Compositional Data Analysis in Practice, Chapman & Hall / CRC. \cr
}
\author{
Michael Greenacre
}
\examples{
data(cups)
apply(cups, 2, sum)
cups <- CLOSE(cups)
apply(cups, 2, sum)
}


