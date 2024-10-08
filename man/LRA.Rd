\name{LRA}
\alias{LRA}
\title{
Logratio analysis
}
\description{
Computation of weighted or unweighted logratio analysis of a samples-by-parts compositional data table.
}
\usage{
LRA(data, nd = 2, weight = TRUE, suprow = NA, row.wt = NA, amalg = NA, supamalg = FALSE)
}
\arguments{
  \item{data}{
A data frame or matrix of compositional data, with no zero values
}
 \item{nd}{
Number of dimensions for summary solution if not 2 (default)
}
  \item{weight}{
TRUE (default) for part weighting, FALSE for unweighted analysis, or a vector of user-defined part weights
}
 \item{suprow}{
Indices of rows that are supplementary points
}
 \item{row.wt}{
Optional user-defined set of positive weights for the rows (samples) (default: equal weights)
}
 \item{amalg}{
Optional list of amalgamated parts 
}
 \item{supamalg}{
FALSE (default) when amalgamations are active and their subparts supplementary, TRUE when amalgamations are supplementary and their parts active  
}
}
\details{
The function \code{LRA} computes a log-ratio analysis of a table of compositional data based on the singular value decomposition. 
By default the weighted log-ratio analysis is computed (Greenacre & Lewi 2009). 
For the unweighted logratio analysis (Aitchison & Greenacre 2002), specify the option \code{weight=FALSE}. 

User-specified weights can be supplied, for the rows and/or the columns.
Usually row weights are not specified, and are equal unless intentional weighting of the samples is desired.
Default column weights (if \code{weight = TRUE}) are the part means of the true compositional table, thus summing to 1.
User-specified part weights can be provided using the same \code{weight} option. 

Supplementary rows can be declared (also known as passive points) -- these do not contribute to the solution but are positioned on the solution axes.

Amalgamations can be defined and can either replace their constituent parts (default) or be declared supplementary using the \code{supamalg} option: \code{supamalg = FALSE} (default), \code{= TRUE} if all declared amalgamations are supplementary.

The function borrows the structure and functions of the \code{ca} package, which is required, and produces a \code{ca} object, and the same \code{print}, \code{summary} and \code{plot} methods can be used, as for a \code{ca} object.}
\value{
  \item{sv        }{Singular values}
  \item{nd        }{Number of dimensions in solution results}
  \item{rownames  }{Row names}
  \item{rowmass   }{Row weights}
  \item{rowdist   }{Row logratio distances to centroid}
  \item{rowinertia}{Row inertias}
  \item{rowcoord  }{Row standard coordinates}
  \item{rowpcoord }{Row principal coordinates}
  \item{rowsup    }{Indices of row supplementary points}
  \item{colnames  }{Column names}
  \item{colmass   }{Column weights}
  \item{coldist   }{Column logratio distances to centroid}
  \item{colinertia}{Column inertias}
  \item{colcoord  }{Column standard coordinates}
  \item{colpcoord }{Column principal coordinates}
  \item{N         }{The compositional data table}
}
\references{
Aitchison, J. and Greenacre, M. (2002), Biplots of compositional data, Applied Statistics 51, 375-392.\cr
Greenacre, M. and Lewi, P.J. (2009), Distributional equivalence and subcompositional coherence in the analysis of compositional data, contingency tables and ratio scale measurements. Journal of Classification 26, 29-54.
Greenacre, M. (2020), Amalgamations are valid in compositional data analysis, can be used in agglomerative clustering, and their logratios have an inverse transformation. Applied Computing and Geosciences 5, 100017.
}
\author{
Michael Greenacre
}
\seealso{
\code{\link[ca]{plot.ca}}, \code{\link[ca]{summary.ca}}, \code{\link[ca]{print.ca}}
}
\examples{
# (weighted) LRA of the RomanCups data set, showing default symmetric map
data("cups")
PLOT.LRA(LRA(cups))
# (unweighted) LRA of the RomanCups data set, showing default symmetric map
# the solution is completely different
PLOT.LRA(LRA(cups, weight=FALSE))
}
\keyword{multivariate}
