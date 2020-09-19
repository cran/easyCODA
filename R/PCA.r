PCA <- function(data, nd = 2, weight = TRUE, row.wt = NA, suprow = NA) {
  # data          interval-scaled data table (e.g. logratios)
  # nd            number of dimensions for summary output, 2 by default
  # weight        = TRUE (default) when weights are in logratio list object
  #               = FALSE for unweighted analysis
  #               = vector of user-specified column weights
  # row.wt        optional weights for rows
  # suprow        columns of supplementary (passive) rows

  # 15/5/2019  fixed bug in supplementary row option, disabled supcol since it makes little sense for CoDa
  # 07/02/2020 fixed bug in row.wt definition

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
  Xact <- X
  if(!is.na(suprow[1])) { Xact <- X[-suprow,] }
  #### option omitted  if(!is.na(supcol[1])) { Xact <- X[,-supcol] }

  if(!weight[1]) cm <- rep(1/ncol(Xact), ncol(Xact))
  if(length(weight) > 1) {
    if(length(weight) != ncol(data)) stop("Number of specified column weights not equal to number of columns of data")
    if(sum(weight <= 0) > 1) stop("Column weights have to be all positive")
    cm <- weight / sum(weight)
  }
  if(length(data.wt) == ncol(data)) cm <- data.wt / sum(data.wt)
  if((as.numeric(weight[1]) == 1) & (!is.list(foo))) stop("weight = TRUE but data is not a logratio list object with weights")

  # row masses are equal unless specified
  # row masses of supplementary points are 0
  rm <- rep(0, nrow(X))
  if(is.na(suprow[1]))  {
    rm <- rep(1/nrow(X), nrow(X))
    rmact <- rm
  }
  if(!is.na(suprow[1])) {
    rm[(1:nrow(X))[-suprow]] <- rep(1/nrow(Xact), nrow(Xact))
    rmact <- rm[-suprow]
  }

  # for specified weights, a list with as many rows as data must be supplied (incl. 0s for suprow)
  if(!is.na(row.wt[1])) {
    if(length(row.wt) != nrow(data)) stop("Number of specified row weights not equal to number of rows of data")
    if(sum(row.wt < 0) > 0) stop("Row weights have to be all nonnegative")
    if(!is.na(suprow[1]) & (sum(row.wt[suprow]!=0)>0)) stop("Supplementary row weights should be 0")
    rm <- rep(0, nrow(X))
    if(is.na(suprow[1]))  rmact <- row.wt / sum(row.wt)
    if(!is.na(suprow[1])) rmact[(1:nrow(X))[-suprow]] <- row.wt / sum(row.wt)
  }
  Y    <- Xact
  mc   <- t(Y) %*% as.vector(rmact)
  Y    <- Y - rep(1,nrow(Xact)) %*% t(mc)
  Z    <- diag(sqrt(rmact)) %*% Y %*% diag(sqrt(cm))
  svdZ <- svd(Z)

  ndmax    <- min(nrow(Xact), ncol(Xact))
  data.rsc <- (diag(1/sqrt(rmact)) %*% svdZ$u)[,1:ndmax]
  data.csc <- (diag(1/sqrt(cm)) %*% svdZ$v)[,1:ndmax]

  # supplementary rows but no supplementary columns (latter option non-existent)
  #  if(!is.na(suprow[1]) & is.na(supcol[1])) {
  if(!is.na(suprow[1])) {

    # supplementary rows; mc is the mean column vector of the active set
    Ysup     <- as.matrix(X[suprow,])
    if(length(suprow)==1) Ysup <- matrix(Ysup - mc, nrow=1)
    if(length(suprow)>1)  Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mc)
    Ysup.rsc <- Ysup %*% diag(sqrt(cm)) %*% svdZ$v[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax])

    # insert supplementary row information in correct places
    foo           <- matrix(0, nrow = nrow(data), ncol = ndmax)
    foo[-suprow,] <- data.rsc
    foo[suprow,]  <- Ysup.rsc
    data.rsc      <- foo
  }

  #### supplementary columns but no supplementary rows
  ###  if(is.na(suprow[1]) & !is.na(supcol[1])) {

  ### supplementary columns, work with matrix transpose
  ###    Ysup     <- t(as.matrix(Xact[,supcol]))
  ###    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mr)
  ###    mc       <- Ysup %*% as.vector(rm)
  ###    Ysup     <- Ysup - mc %*% t(rep(1, ncol(Ysup)))
  ###    Ysup.csc <- Ysup %*% diag(sqrt(rm)) %*% svdZ$u[,1:ndmax] %*% diag(1/svdZ$d[1:ndmax])

  #### insert supplementary row information in correct places
  ###    foo           <- matrix(0, nrow = ncol(data), ncol = ndmax)
  ###    foo[-supcol,] <- data.csc
  ###    foo[supcol,]  <- Ysup.csc
  ###    data.csc      <- foo
  ###  }

  # principal coordinates
  data.rpc <- data.rsc %*% diag(svdZ$d[1:ndmax])
  data.cpc <- data.csc %*% diag(svdZ$d[1:ndmax])

  # substitute results in ca object - first make dummy run
  data.pca             <- ca(abs(X+1), suprow=suprow)
  data.pca$nd          <- nd
  data.pca$rowsup      <- suprow
  data.pca$colsup      <- NA
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
  data.pca$call        <- match.call()
  return(data.pca)
}