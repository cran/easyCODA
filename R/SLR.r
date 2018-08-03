SLR <- function(data, numer=NA, denom=NA, weight=TRUE) {
# computes an amalgamation logratio (ILR) --- "S" stands for "summed"
# numer = set of parts in numerator, denom = set of parts in denominator
# weight = FALSE (unweighted), = TRUE (weighted by column means)
#        = (vector of prespecified weights)
  if(sum(data==0) > 0) stop("Error: some data values are zero")
  data <- as.matrix(data / apply(data, 1, sum))
  if(!weight[1]) weights <- rep(1/ncol(data), ncol(data))
  if(weight[1])  weights <- apply(data, 2, mean)
  if(length(weight) == ncol(data)) {
    if(sum(weight<=0) > 0) stop("Error: some weights zero or negative")
    if(sum(weight)!=1) print("Sum of weights not exactly 1, but are rescaled")
    weights <- weight / sum(weight)
  }
  if(is.na(numer[1]) | is.na(denom[1])) stop("Numerator and/or denominator not specified")
  if(length(intersect(denom,numer)) > 0) stop("Error: numerator and denominator intersect")
  if(length(numer) == 1) num <- data[,numer]
  if(length(numer) > 1) num <- apply(data[,numer], 1, sum)
  if(length(denom) == 1) den <- data[,denom]
  if(length(denom) > 1) den <- apply(data[,denom], 1, sum)
  num.wt <- sum(weights[numer])
  den.wt <- sum(weights[denom])
# coefficient not necessary for amalgamation logratios...
# amalglr <- sqrt(num.wt * den.wt / (num.wt + den.wt)) * log(num/den)
  amalglr <- log(num/den)
  amalglr.weight <- num.wt * den.wt
  names(amalglr.weight) <- paste(paste(colnames(data)[numer], collapse="&"), paste(colnames(data)[denom], collapse="&"), sep="/")
  return(list(LR=amalglr, LR.wt=amalglr.weight))
}
