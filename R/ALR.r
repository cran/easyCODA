ALR <- function(data, denom=ncol(data), weight=TRUE, stats=FALSE) {
# Computes additive logratios (ALRs)
# denom  = reference part (default, last part)
# weight = TRUE (weighted by column means, default),
#        = FALSE (unweighted), 
#        = (vector of prespecified positive weights)
# stats  = TRUE (compute means, variances and total variance)
#        = FALSE (no statistics output) 
  if(sum(data==0) > 0) stop("Error: some data values are zero")
  data <- as.matrix(data / apply(data, 1, sum))
  if(!weight[1]) weights <- rep(1/ncol(data), ncol(data))
  if(weight[1])  weights <- apply(data, 2, mean)
  if(length(weight) == ncol(data)) {
    if(sum(weight<=0) > 0) stop("Error: some weights zero or negative")
    if(sum(weight)!=1) print("Sum of weights not exactly 1, but are rescaled")
    weights <- weight / sum(weight)
  }
  alr <- log(data[,-denom]/data[,denom])
  colnames(alr) <- paste(colnames(data)[-denom], colnames(data)[denom], sep="/")
  alr.weights <- weights[-denom] * weights[denom]
  names(alr.weights) <- colnames(alr)
  if(!stats) return(list(LR=alr, LR.wt=alr.weights, denom=denom, part.names=colnames(data), part.wt=weights))
  if(stats) {
    means <- apply(alr, 2, mean)
    vars  <- apply(alr, 2, var)
    totvar <- sum(vars)
    return(list(LR=alr, LR.wt=alr.weights, denom=denom, part.names=colnames(data), part.wt=weights,
                means=means, vars=vars, totvar=totvar))
  }
}

