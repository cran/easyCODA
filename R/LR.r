LR <- function(data, ordering=1:ncol(data), weight=TRUE) {
# computes all logratios (LRs)
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
  lr <- matrix(0, nrow(data), ncol(data) * (ncol(data)-1) / 2)
  data <- data[, ordering]
  colnames(lr) <- 1:(ncol(data) * (ncol(data)-1) / 2)
  lr.weights <- rep(0, ncol(data) * (ncol(data)-1) / 2)
  k <- 1
  for(j in 1:(ncol(data)-1)) {
    for(jj in (j+1):ncol(data)) {
      lr[,k] <- log(data[,j]/data[,jj])
      colnames(lr)[k] <- paste(colnames(data)[j], colnames(data)[jj], sep="/")
      lr.weights[k] <- weights[j] * weights[jj]   
      names(lr.weights)[k] <- colnames(lr)[k]
      k <- k + 1
    }
  }
  rownames(lr) <- rownames(data)
  return(list(LR=lr, LR.wt=lr.weights))
}


