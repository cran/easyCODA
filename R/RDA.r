RDA <- function(data, cov=NA, nd = NA, weight = TRUE, suprow = NA,
                row.wt = NA) {
# data          interval-scaled data table (e.g. logratios)
# cov           list of covariates for constraining solution
# nd            number of dimensions for summary output, by default the number of constrained dimensions
# weight        = TRUE (default) when weights are in logratio list object
#               = FALSE for unweighted analysis
#               = vector of user-specified column weights   
# suprow        columns of supplementary (passive) rows
# row.wt        optional weights for rows


  supcol <- NA   # this option not active yet
  data.LR <- data
  data.wt <- weight   # either TRUE/FALSE or equal to list of weights
  cm <- 0

# first check if there's anything in data$LR.wt, otherwise regular matrix
  if(!is.matrix(data)) {
    if(length(data$LR.wt) > 0) {
      data.wt <- data$LR.wt
      data.LR <- data$LR
      if(weight[1]) cm <- data.wt
    }
  }

# case weight=TRUE but weights were not in data object
  if((weight[1]) & (length(data.wt)==1)) stop("weight = TRUE but data is not a list object with weights")

# set up matrix X for PCA
  X    <- as.matrix(data.LR)
  Xsup <- X
  if(!is.na(suprow[1])) { Xsup <- X[-suprow,] }
  if(!is.na(supcol[1])) { Xsup <- X[,-supcol] }
  if(!is.na(supcol[1]) & length(cm)>0) cm <- cm[-supcol]

# case weight=FALSE
  if(!weight[1]) cm <- rep(1/ncol(Xsup), ncol(Xsup))

# case weight=column weights
  if(length(weight) > 1) {
    if(length(weight) != ncol(data.LR)) stop("Number of specified column weights not equal to number of columns of data")
    if(sum(weight <= 0) > 1) stop("Column weights have to be all positive")
    cm <- weight
  }

  rm <- rep(1/nrow(Xsup), nrow(Xsup))
  if(!is.na(row.wt[1])) {
    if(length(row.wt) != nrow(data)) stop("Number of specified row weights not equal to number of rows of data")
    if(sum(row.wt <= 0) > 1) stop("Row weights have to be all positive")
    rm <- row.wt / sum(row.wt)
  }

# RDA computations

# response matrix, centred, weighted
  Y    <- Xsup
  mc   <- t(Y) %*% as.vector(rm)
  Y    <- Y - rep(1,nrow(Xsup)) %*% t(mc)
  S    <- diag(sqrt(rm)) %*% Y %*% diag(sqrt(cm))

# prepare predictor matrix, which can be a mixture of continuous, crisp and fuzzy categorical variables
# index values of the continuous and categorical (dummy & fuzzy variables)  

  Z    <- as.matrix(cov)
  Dr   <- diag(rm)
  one  <- as.matrix(rep(1,nrow(Z)))
  Z    <- Z-one%*%t(rm)%*%Z
  wvar <- diag(as.matrix(t(Z)%*%Dr%*%Z))
  Z    <- Z%*%sqrt(diag(1/wvar))

# weighted correlation matrix
  R    <- as.matrix(t(Z)%*%Dr%*%Z)

# inverse of weighted correlation matrix
  R.svd <- svd(R)
  k     <- length(which(R.svd$d > 1.e-10))
  invR  <- as.matrix(R.svd$u[,1:k]) %*% diag(1/R.svd$d[1:k], nrow=k, ncol=k) %*% t(as.matrix(R.svd$v[,1:k]))

# projection matrix and projected responses
  Proj <- sqrt(Dr) %*% Z %*% invR %*% t(Z) %*% sqrt(Dr)
  A    <- Proj %*% S
  svdA <- svd(A)

# standard and principal coordinates
  ndmax    <- sum(svdA$d > 1.e-10)
  data.rsc <- (diag(1/sqrt(rm)) %*% as.matrix(svdA$u)[,1:ndmax])
  data.csc <- (diag(1/sqrt(cm)) %*% as.matrix(svdA$v)[,1:ndmax])
  data.rpc <- data.rsc %*% diag(svdA$d[1:ndmax], nrow=ndmax, ncol=ndmax)
  data.cpc <- data.csc %*% diag(svdA$d[1:ndmax], nrow=ndmax, ncol=ndmax)

# supplementary rows but no supplementary columns
  if(!is.na(suprow[1]) & is.na(supcol[1])) { 

# supplementary rows
    Ysup     <- as.matrix(Xsup[suprow,])
    Ysup     <- Ysup - rep(1, nrow(Ysup)) %*% t(mc)
    mr       <- Ysup %*% as.vector(cm)
    Ysup     <- Ysup - mr %*% t(rep(1, ncol(Ysup)))
    Ysup.rsc <- Ysup %*% diag(sqrt(cm)) %*% svdA$v[,1:ndmax] %*% diag(1/svdA$d[1:ndmax], nrow=ndmax, ncol=ndmax)

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
    Ysup.csc <- Ysup %*% diag(sqrt(rm)) %*% svdA$u[,1:ndmax] %*% diag(1/svdA$d[1:ndmax], nrow=ndmax, ncol=ndmax)

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

  
# coordinates of covariates
  cov.reg <- matrix(0, nrow=ncol(cov), ncol=ndmax)
  for(j in 1:ncol(cov)) cov.reg[j,] <- lm(scale(cov[,j]) ~ data.rsc)$coefficients[2:(ndmax+1)]

# substitute results in ca object - first make dummy run
  data.rda             <- ca(abs(Xsup))  
  data.rda$nd          <- ndmax
  data.rda$rowsup      <- suprow
  data.rda$colsup      <- supcol
  data.rda$sv          <- svdA$d[1:ndmax]
  data.rda$rowmass     <- rm
  data.rda$colmass     <- cm
  data.rda$rowcoord    <- data.rsc
  data.rda$colcoord    <- data.csc
  data.rda$rowpcoord   <- data.rpc
  data.rda$colpcoord   <- data.cpc
  data.rda$covcoord    <- cov.reg
  data.rda$covnames    <- colnames(cov)
  data.rda$rowdist     <- sqrt(apply(data.rpc^2, 1, sum))
  data.rda$coldist     <- sqrt(apply(data.cpc^2, 1, sum))
  data.rda$rowinertia  <- apply(diag(rm) %*% as.matrix(data.rpc^2), 1, sum)
  data.rda$colinertia  <- apply(diag(cm) %*% as.matrix(data.cpc^2), 1, sum)
  data.rda$N           <- data.LR
  data.rda$cov         <- cov
  return(data.rda)
  }
