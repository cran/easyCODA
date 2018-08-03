\name{DUMMY}
\alias{DUMMY}
\title{
Dummy variable (indicator) coding
}
\description{
Convert categorical variable to dummy (0/1) coding
}
\usage{
DUMMY(x, catnames=NA)
}
\arguments{
  \item{x}{
Variable (vector) of categorical data to be coded
}
 \item{catnames}{
Category names
}
}
\details{
The function \code{DUMMY} takes a categorical variable and converts it to a set of dummy variables (zeros and ones), where the ones indicate the corresponding category. There are as many columns in the result as there are unique categories in the input vector. 
}
\author{
Michael Greenacre
}
\examples{
# Indicator (dummy) coding of sex in FishMorphology data set
data(fish)
sex   <- fish[,1]
sex.Z <- DUMMY(sex, catnames=c("F","M")) 
}
\keyword{coding}
