PLR <- function(data, ordering=1:ncol(data), weight=TRUE) {
# computes the pivot logratios (PLRs)
# ordering is a permutation of the columns, by default the original ordering
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
  if(sum(intersect(1:ncol(data), ordering) == 1:ncol(data)) != ncol(data)) 
     stop("ordering is not a permutation of the columns of data matrix")
  plr <- matrix(0, nrow(data), ncol(data)-1)
  data <- data[, ordering]
  colnames(plr) <- 1:(ncol(data)-1)
  plr.weights <- rep(0, ncol(data)-1)
  for(j in 1:(ncol(data)-1)) {
    ilr <- ILR(data, numer=j, denom=(j+1):ncol(data))
    plr[,j] <- ilr$LR
    colnames(plr)[j] <- names(ilr$LR.wt)
    plr.weights[j] <- ilr$LR.wt    
    names(plr.weights)[j] <- names(ilr$LR.wt)
  }
  rownames(plr) <- rownames(data)
  return(list(LR=plr, LR.wt=plr.weights))
}
