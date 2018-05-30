LR.VAR <- function(LRdata, row.wt = NA, weight=TRUE, vars=FALSE) {

# LRdata = matrix of logratios, can also be a single vector of logratios
#   should be output from one of the functions ALR, CLR, ILR, SLR, PLR or LR,
#   which include weights in the objective;
# alternatively a set of any data, where column weights are specified explicitly
#   in the option weight
# row.wt = row weights, summing to 1 
# weight = TRUE (weighted by LR weights, default) 
#        = FALSE (unweighted)
#        = (vector of prespecified weights)
# vars   = TRUE (include in a list the individual variances) 
#        = FALSE (don't include all individual variances, default)
#        

  if(!is.list(LRdata)) {
    LRfoo <- as.matrix(LRdata)
    weights <- 1
  }
  if(is.list(LRdata)) {
    LRfoo <- LRdata$LR
    if(!weight[1]) weights <- rep(1/ncol(LRfoo), ncol(LRfoo))
    if(weight[1] & (ncol(LRfoo)>1)) weights <- LRdata$LR.wt
    if(length(weight) == ncol(LRfoo)) {
      if(sum(weight<=0) > 0) stop("Error: some weights zero or negative")
      if(sum(weight)!=1) print("Sum of column weights not exactly 1, but are rescaled")
      weights <- weight / sum(weight)
    }
  }
  if(is.na(row.wt[1])) row.wt <- rep(1/nrow(LRfoo), nrow(LRfoo))
  if(length(row.wt) != nrow(LRfoo)) stop("Error: row weights not the same number as rows of data")
  if(sum(row.wt)!=1) print("Sum of row weights not exactly 1, but are rescaled")
  row.wt <- row.wt/sum(row.wt)
# only one variable, single variance
  if(ncol(LRfoo) == 1) {
    LRmean    <- sum(LRfoo * row.wt)
    LRtotvar  <- sum((LRfoo-LRmean)^2 * row.wt)
  }
  if(ncol(LRfoo) > 1) {
# column-centre the data matrix, many variances
    LRfoo <- sweep(LRfoo, 2, apply(LRfoo * row.wt, 2, sum)) 
    LRtotvar <- sum(diag(row.wt) %*% LRfoo^2 %*% diag(weights))
    if(vars) {
      LRvars <- apply(diag(row.wt) %*% LRfoo^2 %*% diag(weights), 2, sum)
      names(LRvars) <- colnames(LRfoo)
    }
  }
  if(!vars) return(LRtotvar)
  if(vars)  return(list(LRtotvar=LRtotvar, LRvars=LRvars))
}

