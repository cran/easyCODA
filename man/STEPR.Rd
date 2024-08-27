\name{STEPR}
\alias{STEPR}
\title{
Stepwise selection of pairwise logratios for generalized linear modelling
}
\description{
Three different algorithms for selecting pairwise logratios that best explain/predict a response variable, which could be continuous, binary or count
}
\usage{
STEPR(data, y, method = NA, family = "gaussian", nsteps = ncol(data)-1, 
      top = 1, previous = NA, criterion = "Bonferroni", alpha = 0.05, 
      previousparts=NA, denom=NA)
}
\arguments{
  \item{data}{
A data frame or matrix of compositional data on which the pairwise logratios will be constructed and selected
}
 \item{y}{
The response variable: a numeric variable for regression (default)), a binary factor for logistic regression or a numeric count for Poisson regression
}
 \item{method}{
The selection method: 1 (unrestricted selection of logratios), 2 (restricted to non-overlapping parts), 3 (additive logratios)
}
 \item{family}{
The distribution used in the generalized linear model family: "gaussian" (default, for multiple regression), "binomial" (for logistic regression of binary response), or "poisson" (for Poisson regression)
}
  \item{nsteps}{
The maximum number of steps taken, by default one less than the number of parts
}
 \item{top}{
When one step is taken (nsteps=1), the ordered list of top logratios with the highest improvements in the likelihood function, for selection based on domain knowledge 
}
 \item{previous}{	
For specifying variable(s) to be included before stepwise selection takes place; these can be non-compositional variables and/or specific pairwise logratios computed in previous runs of STEPR or by hand; the matrix (or vector for a single variable) of values must be supplied
}
 \item{criterion}{
Criterion for stopping the stepwise selection: "Bonferroni" (default), "AIC", "BIC", or NA for no stopping until maximum specified or permissible logratios entered
}
 \item{alpha}{	
Overall significance level (default is 0.05)
}
 \item{previousparts}{
(For method 2) The sequence numbers of the logratios, if any, forced in using the previous option
}
 \item{denom}{
(For method 3) The sequence number of the part used in denominator; for use when additive logratios are forced in using previous option or to select a set of additive logratios with specific reference from the start
}
}
\details{
The function \code{STEPR} performs stepwise selection of pairwise logratios with the objective of explaining/predicting a response variable, in the framework of generalized linear modelling where the response can be numeric continuous (regression analysis), or a binary factor (logistic regression), or a numeric count (Poisson regression). The corresponding family option has to be indicated if the regression is logistic or Poisson. The different method options for the stepwise selection are \code{method = 1} (unrestricted selection of logratios, any logratios can be selected irrespective of the previous ones), \code{method = 2} (restricted to non-overlapping parts, each part participates at most in one logratio, so that parts in previously selected logratios are excluded in subsequent steps; logratio effects can be interpreted as under orthogonality), \code{method = 3} (additive logratios, only logratios with the same denominator as the first selected logratio are c onsidered; the result is an additive logratio transformation on a subcomposition) Three alternative stopping criteria can be specified, otherwise the procedure executes as many steps as the value of nsteps. These are (in increasing strictness), "AIC", "BIC" and "Bonferroni" (the default).
}
\value{
  \item{rationames        }{Names of the selected logratios}
  \item{ratios       }{The sequence numbers of the selected parts in each ratio}
  \item{logratios    }{Matrix of selected logratios}
  \item{logLik        }{The -2*log-likelihood sequence for the steps}
  \item{deviance      }{The deviance sequence for the steps}
  \item{AIC    }{The AIC sequence for the steps}
  \item{BIC   }{The BIC sequence for the steps}
  \item{Bonferroni}{The Bonferroni sequence for the steps}
  \item{null.deviance       }{The null deviance for the regression}
  
(Notice that for \code{logLik}, \code{AIC}, \code{BIC} and \code{Bonferroni}, the values for one more step are given, so that the stopping point can be confirmed.)  \cr
  
And the following if \code{top > 1}: \cr

  \item{ratios.top}{The top ratios and the sequence numbers of their parts}
  \item{logratios.top}{The matrix of top logratios}
  \item{logLik.top}{The set of top -2*log-likelihoods}
  \item{deviance.top}{The set of top deviances}
  \item{AIC.top}{The set of top AICs}
  \item{BIC.top}{The set of top BICs}
  \item{Bonferroni.top}{The set of top Bonferronis}
}
\references{
Coenders, G. and Greenacre, M. (2021), Three approaches to supervised learning for compositional data with pairwise logratios. aRxiv preprint. URL:https://arxiv.org/abs/2111.08953 \cr
Coenders, G. and Pawlowsky-Glahn, V. (2020), On interpretations of tests and effect sizes in regression models with a compositional predictor. SORT, 44:201-220 \cr
Greenacre, M. (2021), Compositional data analysis, Annual Review of Statistics and its Application, 8: 271-299
}
\author{
Michael Greenacre
}
\seealso{
\code{\link{ALR}}, \code{\link{STEP}}, \code{\link{glm}}
}
\examples{
# For the fish morphometric data, first close (normalize, although not necessary) 
# then loop over the 26*25/2 = 325 possible logratios stepwise
data(fish)
habitat <- fish[,2]
morph <- CLOSE(fish[,4:29])
# predict habitat binary classification from morphometric ratios
fish.step1 <- STEPR(morph, as.factor(habitat), method=1, family="binomial")
# [1] "Criterion increases when 3-th ratio enters"
fish.step1$names
# [1] "Bac/Hg" "Hw/Jl"
# perform logistic regression with selected logratios
fish.glm   <- glm(as.factor(habitat) ~ fish.step1$logratios, family="binomial")
summary(fish.glm)
fish.pred1  <- predict(fish.glm)
table(fish.pred1>0.5, habitat)
#     habitat
#        1  2
# FALSE 56 11
# TRUE   3  5
# (Thus 61/75 correct predictions)
#
# force the sex variable in at the first step before selecting logratios
# and using more strict Bonferroni default
sex <- as.factor(fish[,1]) 
fish.step2 <- STEPR(morph, as.factor(habitat), method=1, previous=sex, family="binomial")
# [1] "Criterion increases when 3-th ratio enters"
fish.step2$names
# [1] "Bac/Hg" "Hw/Jl"
# perform logistic regression with sex and selected logratios
fish.glm   <- glm(as.factor(habitat) ~ sex + fish.step2$logratios, family="binomial")
summary(fish.glm)
# (sex not significant)
#
# check the top 10 ratios at Step 1 to allow domain knowledge to operate
fish.step3 <- STEPR(morph, as.factor(habitat), method=1, nsteps=1, top=10, family="binomial")
cbind(fish.step3$ratios.top, fish.step3$BIC.top)
#        row col         
# Bac/Hg   8  19 67.93744
# Bp/Hg    7  19 69.87134
# Jl/Hg    6  19 70.31554
# Jw/Bp    5   7 71.53671
# Jw/Jl    5   6 71.57122
# Jw/Bac   5   8 71.69294
# Fc/Hg   10  19 72.38560
# Hw/Bac   1   8 73.25325
# Jw/Fc    5  10 73.48882
# Hw/Bp    1   7 73.55621
# Suppose 5th in list, Jw/Jl (Jaw width/Jaw length), preferred at the first step
fish.step4 <- STEPR(morph, as.factor(habitat), method=1, 
                    previous=fish.step3$logratios.top[,5], family="binomial")
# [1] "Criterion increases when 2-th ratio enters"
fish.step4$names
# [1] "Bac/Hg"
# So after Jw/Jl forced in only Bac/Hg enters, the best one originally
}
\keyword{multivariate}
