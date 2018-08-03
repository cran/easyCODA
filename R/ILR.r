ILR <- function(data, numer=NA, denom=NA, weight=TRUE) {
# computes one isometric logratio (ILR)
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
  if(length(numer) == 1) num <- log(data[,numer]) * weights[numer]
  if(length(numer) > 1) num <- apply(log(data[,numer]) %*% diag(weights[numer]), 1, sum)
  num <- num / sum(weights[numer])
  if(length(denom) == 1) den <- log(data[,denom]) * weights[denom]
  if(length(denom) > 1) den <- apply(log(data[,denom]) %*% diag(weights[denom]), 1, sum)
  den <- den / sum(weights[denom])
  num.wt <- sum(weights[numer])
  den.wt <- sum(weights[denom])
  ilr <- sqrt(num.wt * den.wt / (num.wt + den.wt)) * (num - den)
  ilr.weight <- num.wt * den.wt
  names(ilr.weight) <- paste(paste(colnames(data)[numer], collapse="&"), paste(colnames(data)[denom], collapse="&"), sep="/")
  return(list(LR=ilr, LR.wt=ilr.weight))
}
