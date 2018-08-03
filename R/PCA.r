PCA <- function(data, nd = 2, weight = TRUE, row.wt = NA, suprow = NA, supcol = NA) {
# data          interval-scaled data table (e.g. logratios)
# nd            number of dimensions for summary output, 2 by default
# weight        = TRUE (default) when weights are in logratio list object
#               = FALSE for unweighted analysis
#               = vector of user-specified column weights  
# row.wt        optional weights for rows 
# suprow        columns of supplementary (passive) rows
# supcol        columns of supplementary (passsive) columns

# first check if data is a dataframe, if so make it a matrix
  if(is.data.frame(data)) data <- as.matrix(data)
   
# then check if data is a list
  foo <- data
  data.wt <- 0
  if(is.list(foo)) {
    data <- foo$LR
    data.wt <- foo$LR.wt
  }

# do PCA
  X    <- as.matrix(data)
  Xsup <- X
  if(!is.na(suprow[1])) { Xsup <- X[-suprow,] }
  if(!is.na(supcol[1])) { Xsup <- X[,-supcol] }

  if(!weight[1]) cm <- rep(1/ncol(Xsup), ncol(Xsup))
  if(length(weight) > 1) {
    if(length(weight) != ncol(data)) stop("Number of specified column weights not equal to number of columns of data")
    if(sum(weight <= 0) > 1) stop("Column weights have to be all positive")
    cm <- weight / sum(weight)
  }
  if(length(data.wt) == ncol(data)) cm <- data.wt / sum(data.wt)
  if((as.numeric(weight[1]) == 1) & (!is.list(foo))) stop("weight = TRUE but data is not a list object with weights")
  rm <- rep(1/nrow(Xsup), nrow(Xsup))
  if(!is.na(row.wt[1])) {
    if(length(row.wt) != nrow(data)) stop("Number of specified row weights not equal to number of rows of data")
    if(sum(row.wt <= 0) > 1) stop("Row weights have to be all positive")
    rm <- row.wt / sum(row.wt)
  }
  Y    <- Xsup
  mc   <- t(Y) %*% as.vector(rm)
  Y    <- Y - rep(1,nrow(Xsup)) %*% t(mc)
  Z    <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))
  svdZ <- svd(Z)

  ndmax    <- min(nrow(Xsup), ncol(Xsup))
  data.rsc <- (diag(1/sqrt(rm)) %*% svdZ$u)[,1:ndmax]
  data.csc <- (diag(1/sqrt(cm)) %*% svdZ$v)[,1:ndmax]

# supplementary rows but no supplementary columns
  if(!is.na(suprow[1]) & is.na(supcol[1])) { 

# supplementary rows
    Ysup     <- as.matrix(Xsup[suprow,])
    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mc)
    mr       <- Ysup %*% as.vector(cm)
    Ysup     <- Ysup - mr %*% t(rep(1, ncol(Ysup)))
    Ysup.rsc <- Ysup %*% diag(sqrt(cm)) %*% svdZ$v[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax]) 

# insert supplementary row information in correct places
    foo           <- matrix(0, nrow = nrow(data), ncol = ndmax) 
    foo[-suprow,] <- data.rsc
    foo[suprow,]  <- Ysup.rsc 
    data.rsc      <- foo
#    foo           <- rep(0, nrow(data))  
#    foo[-suprow]  <- rm
#    foo[suprow]   <- apply(P[suprow,], 1, sum) / sum(P)
    rm            <- rm
  }

# supplementary columns but no supplementary rows
  if(is.na(suprow[1]) & !is.na(supcol[1])) { 

# supplementary columns, work with matrix transpose
    Ysup     <- t(as.matrix(Xsup[,supcol]))
    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mr)
    mc       <- Ysup %*% as.vector(rm)
    Ysup     <- Ysup - mc %*% t(rep(1, ncol(Ysup)))
    Ysup.csc <- Ysup %*% diag(sqrt(rm)) %*% svdZ$u[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax]) 

# insert supplementary row information in correct places
    foo           <- matrix(0, nrow = ncol(data), ncol = ndmax) 
    foo[-supcol,] <- data.csc
    foo[supcol,]  <- Ysup.csc 
    data.csc      <- foo
#    foo           <- rep(0, ncol(data))
#    foo[-supcol]  <- cm
#    foo[supcol]   <- apply(P[,supcol], 2, sum) / sum(P)
    cm            <- cm
  }

# principal coordinates
  data.rpc <- data.rsc %*% diag(svdZ$d[1:ndmax])
  data.cpc <- data.csc %*% diag(svdZ$d[1:ndmax])
  
# substitute results in ca object - first make dummy run
  data.pca             <- ca(abs(Xsup))  
  data.pca$nd          <- nd
  data.pca$rowsup      <- suprow
  data.pca$colsup      <- supcol
  data.pca$sv          <- svdZ$d[1:ndmax]
  data.pca$rowmass     <- rm
  data.pca$colmass     <- cm
  data.pca$rowcoord    <- data.rsc
  data.pca$colcoord    <- data.csc
  data.pca$rowpcoord   <- data.rpc
  data.pca$colpcoord   <- data.cpc
  data.pca$rowdist     <- sqrt(apply(data.rpc^2, 1, sum))
  data.pca$coldist     <- sqrt(apply(data.cpc^2, 1, sum))
  data.pca$rowinertia  <- apply(diag(rm) %*% (data.rpc^2), 1, sum)
  data.pca$colinertia  <- apply(diag(cm) %*% (data.cpc^2), 1, sum)
  data.pca$N           <- data
  return(data.pca)
  }

