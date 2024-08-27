\name{easyCODA-package}
\alias{easyCODA-package}
\docType{package}
\title{
\packageTitle{easyCODA}
}
\description{
\packageDescription{easyCODA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{easyCODA}
\packageIndices{easyCODA}
}
\author{
\packageAuthor{easyCODA}

Maintainer: \packageMaintainer{easyCODA}
}
\references{
Greenacre, Michael (2018) Compositional Data Analysis in Practice. Chapman & Hall / CRC Press
}
\keyword{package}
\seealso{
\code{\link[ca]{ca}}
}
\examples{
# Roman cups glass compositions
data("cups")
# unweighted logratio analysis
cups.ulra <- LRA(cups, weight=FALSE)
PLOT.LRA(cups.ulra)
# weighted logratio analysis
cups.wlra <- LRA(cups)
PLOT.LRA(cups.wlra)
}
