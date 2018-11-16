LRA <- function(data, nd = 2, weight = TRUE, suprow = NA, row.wt = NA, 
                amalg = NA, supamalg = FALSE) {
# updated 1/9/2018 to fix bug in definition of cm

# data          compositional data table
# nd            number of dimensions for summary output, 2 by default
# weight        = TRUE (default) when weights are column means
#               = FALSE for unweighted analysis
#               = vector of user-specified column weights            
# suprow        columns of supplementary (passive) rows
# row.wt        optional weights for rows
# amalg         optional list of amalgamated parts
# supamalg      = FALSE (default) when amalgamations are active and their subparts supplementary
#               = TRUE when amalgamations are desired to be supplementary and their subparts active

# first sort out amalgamations, if any
  if(!is.na(amalg[1])) {
    foo <- data
    for(a in 1:length(amalg)) {
      if(length(amalg[[a]]) == 1) stop("Amalgamations must have more than 1 part")
      foo <- cbind(foo, apply(foo[,amalg[[a]]], 1, sum))
    }
    colnames(foo)[(ncol(foo)-length(amalg)+1): ncol(foo)] <- names(amalg)
    supcol <- amalg[[1]]
    if(length(amalg) > 1) {
      for(a in 2:length(amalg)) supcol <- c(supcol, amalg[[a]])
    }
    data <- foo
    if(supamalg) supcol <- (ncol(foo)-length(amalg)+1): ncol(foo)
    data.lra <- ca(data, suprow = suprow, supcol = supcol)
  }
  if(is.na(amalg[1])) {
    data.lra <- ca(data, suprow = suprow)
    supcol <- NA
  }
 
# then do actual LRA
  if(sum(data == 0) > 0) {
    stop("There are data zeros in the matrix -- replace with positive values")
  }
  if(sum(data < 0) > 0) {
    stop("There are negative values in the matrix -- this is not allowed for logratio analysis")
  }
  data <- as.matrix(data)
  P    <- data/sum(data)
  Psup <- P
  if(!is.na(suprow[1])) { Psup <- P[-suprow,] ; Psup <- Psup / sum(Psup) }
  if(!is.na(supcol[1])) { Psup <- P[,-supcol] ; Psup <- Psup / sum(Psup) }
  rowsums     <- apply(Psup, 1, sum)
  rowsums.dev <- diff(range(rowsums)) / min(rowsums)
  if(rowsums.dev > 0.01) {
    message("Warning: The row sums of this matrix are not constant, so this may not be a compositional data matrix")
    message("         Nevertheless, the method is still valid and continues since all the data values are positive")
    message("         If this was not intended, please supply a valid compositional table with constant row sums")
  }
  rm   <- apply(Psup,1,sum)
  cm   <- apply(Psup,2,sum)
  if(!weight[1]) cm <- rep(1/ncol(Psup), ncol(Psup))
  if(length(weight) > 1) {
    if(length(weight) != ncol(data)) stop("Number of specified column weights not equal to number of columns of data")
    if(sum(weight <= 0) > 1) stop("Column weights have to be all positive")
    cm <- weight / sum(weight)
  }
  if(!is.na(row.wt[1])) {
    if(length(row.wt) != nrow(data)) stop("Number of specified row weights not equal to number of rows of data")
    if(sum(row.wt <= 0) > 1) stop("Row weights have to be all positive")
    rm <- row.wt / sum(row.wt)
  }
  Y    <- as.matrix(log(Psup))
  mc   <- t(Y) %*% as.vector(rm)
  Y    <- Y - rep(1,nrow(Psup)) %*% t(mc)
  mr   <- Y %*% as.vector(cm)
  Y    <- Y - mr %*% t(rep(1,ncol(Psup)))
  Z    <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))
  svdZ <- svd(Z)

  ndmax    <- min(nrow(Psup), ncol(Psup)) - 1
  data.rsc <- (diag(1/sqrt(rm)) %*% svdZ$u)[,1:ndmax]
  data.csc <- (diag(1/sqrt(cm)) %*% svdZ$v)[,1:ndmax]

# supplementary rows but no supplementary columns
  if(!is.na(suprow[1]) & is.na(supcol[1])) { 

# supplementary rows
    Ysup     <- as.matrix(log(P[suprow,]))
    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mc)
    mr       <- Ysup %*% as.vector(cm)
    Ysup     <- Ysup - mr %*% t(rep(1, ncol(Ysup)))
    Ysup.rsc <- Ysup %*% diag(sqrt(cm)) %*% svdZ$v[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax]) 

# insert supplementary row information in correct places
    foo           <- matrix(0, nrow = nrow(data), ncol = ndmax) 
    foo[-suprow,] <- data.rsc
    foo[suprow,]  <- Ysup.rsc 
    data.rsc      <- foo
    foo           <- rep(0, nrow(data))
    foo[-suprow]  <- rm
    foo[suprow]   <- apply(P[suprow,], 1, sum) / sum(P)
    rm            <- foo
  }

# supplementary columns but no supplementary rows
  if(is.na(suprow[1]) & !is.na(supcol[1])) { 

# supplementary columns, work with matrix transpose
    Ysup     <- t(as.matrix(log(P[,supcol])))
    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mr)
    mc       <- Ysup %*% as.vector(rm)
    Ysup     <- Ysup - mc %*% t(rep(1, ncol(Ysup)))
    Ysup.csc <- Ysup %*% diag(sqrt(rm)) %*% svdZ$u[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax]) 

# insert supplementary column information in correct places
    foo           <- matrix(0, nrow = ncol(data), ncol = ndmax) 
    foo[-supcol,] <- data.csc
    foo[supcol,]  <- Ysup.csc 
    data.csc      <- foo
    foo           <- rep(0, ncol(data))
    foo[-supcol]  <- cm
    foo[supcol]   <- apply(data[,supcol], 2, sum) / sum(data[,-((ncol(data)-length(amalg)+1): ncol(data))])
    cm            <- foo
  }

# principal coordinates (returned)
  data.rpc <- data.rsc %*% diag(svdZ$d[1:ndmax])
  data.cpc <- data.csc %*% diag(svdZ$d[1:ndmax])
  
# substitute results in ca object
  data.lra$nd          <- nd
  data.lra$rowsup      <- suprow
  data.lra$colsup      <- supcol
  data.lra$sv          <- svdZ$d[1:ndmax]
  data.lra$rowmass     <- rm
  data.lra$colmass     <- cm
  data.lra$rowcoord    <- data.rsc
  data.lra$colcoord    <- data.csc
  data.lra$rowpcoord   <- data.rpc
  data.lra$colpcoord   <- data.cpc
  data.lra$rowdist     <- sqrt(apply(data.rpc^2, 1, sum))
  data.lra$coldist     <- sqrt(apply(data.cpc^2, 1, sum))
  data.lra$rowinertia  <- apply(diag(rm) %*% (data.rpc^2), 1, sum)
  data.lra$colinertia  <- apply(diag(cm) %*% (data.cpc^2), 1, sum)
  data.lra$N           <- data
  return(data.lra)
  }
