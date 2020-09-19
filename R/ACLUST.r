ACLUST <- function(data, weight=TRUE, close=TRUE) {
### amalgamation clustering of columns of a compositional data matrix
### default is with weights = part averages summing to 1

### easyCODA preparation 

#  data should be a matrix of compositional data, closed here if close=TRUE
   data <- as.matrix(data)
#  data cannot have any zeros
   if(sum(data==0) > 0) stop("Error: some data values are zero")
   data.pro <- data
   if(close) data.pro <- CLOSE(data)
#  cm are the column massses
   if(length(weight==1) & !weight[1]) cm <- rep(1/ncol(data), ncol(data))
   if(length(weight)==1 & weight[1])  cm <- colMeans(CLOSE(data))

   if(length(weight) > 1) {
     if(length(weight) != ncol(data.pro)) stop("Number of specified column weights not equal to number of columns of data")
     if(sum(weight <= 0) > 0) stop("Column weights have to be all positive")
     cm <- weight / sum(weight)
  }

   data.clr <- CLR(CLOSE(data.pro), weight=cm)$LR
   TOTVAR <- LR.VAR(CLR(CLOSE(data.pro), weight=cm))
#   if(length(weight)==1 & weight[1])  data.clr <- CLR(CLOSE(data.pro), weight=TRUE)$LR
#   if(length(weight) > 1) data.clr <- CLR(CLOSE(data.pro), weight=cm)$LR

   n <- nrow(data.pro)
   m <- ncol(data.pro)

#  data.pro2 is updated in each iteration
   data.pro2 <- data.pro

### diss replaced by easyCODA computation of minimum variance matrix
### (for pairs of *columns*)

   diss <- matrix(99, m, m)
   rownames(diss) <- colnames(diss) <- colnames(data.pro)

   for (j1 in 2:m) {
     for (j2 in 1:(j1-1)) {

# amalgamate the parts in data.pro (which are in columns)
       data.foo      <- data.pro
       foo.j1j2      <- as.numeric(data.pro[,j1] + data.pro[,j2])
       data.foo[,j1] <- foo.j1j2
       data.foo      <- data.foo[,-j2]
#       if(length(weight==1) & !weight[1]) foo.rda <- rda(data.clr ~ CLR(data.foo, weight=FALSE)$LR)
#       if(length(weight)==1 & weight[1])  foo.rda <- rda(data.clr ~ CLR(data.foo, weight=TRUE)$LR)
#       if(length(weight) > 1)             foo.rda       <- rda(data.clr ~ CLR(data.foo, weight=cm)$LR)
#       foo.rda       <- RDA(data.clr, cov=CLR(data.foo, weight=FALSE)$LR, weight=FALSE)
#       diss[j1,j2]   <- diss[j2,j1] <- sum(foo.rda$sv^2)
       if(length(weight==1) & !weight[1]) foo.RDA <- RDA(data.clr, CLR(data.foo, weight=FALSE)$LR, weight=cm)
       if(length(weight==1) & weight[1])  foo.RDA <- RDA(data.clr, CLR(data.foo, weight=TRUE)$LR, weight=cm)
       if(length(weight) > 1)             foo.RDA <- RDA(data.clr, CLR(data.foo, weight=TRUE)$LR, weight=cm)
#       diss[j1,j2]   <- diss[j2,j1] <- (foo.rda$CA$tot.chi/m)*(n-1)/n
       diss[j1,j2]   <- diss[j2,j1] <- TOTVAR - sum(foo.RDA$sv^2)
       }
   }

   flag <- rep(1, m)                          # active/dead indicator
   a <- rep(0, m-1)                           # left subnode on clustering
   b <- rep(0, m-1)                           # right subnode on clustering
   ia <- rep(0, m-1)                          # R-compatible version of a
   ib <- rep(0, m-1)                          # R-compatible version of b
   lev <- rep(0, m-1)                         # level or criterion values
   card <- rep(1, m)                          # cardinalities
   order <- rep(0, m)                         # R-compatible order for plotting
   var.decomp <- rep(0, m-1)                  # variance decomposition
#   nnsnnsdiss <- getnns(diss, flag)           # call to function getnns

### replaced with code to find smallest distance, no error checking needed
### flag all 1s, no need to check
#   steps <- 1

   mm <- rep(0, nrow(diss))
   mmdiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   for (j1 in 1:nrow(diss)) {
          minobs <- -1
          mindis <- MAXVAL
          for (j2 in 1:ncol(diss)) {
              if ( (diss[j1,j2] < mindis) && (j1 != j2) ) {
                 mindis <- diss[j1,j2]
                 minobs <- j2
              }
          }
          mm[j1] <- minobs
          mmdiss[j1] <- mindis
   }

### compatibility with hierclust after "function" call

   mmsmmsdiss <- list(mm=mm, mmdiss=mmdiss)

   clusmat <- matrix(0, m, m)                 # cluster memberships
   for (j in 1:m) clusmat[j,m] <- j           # init. trivial partition

   MAXVAL <- 1.0e12

   for (ncl in (m-1):1) {                      # main loop
       # check for agglomerable pair
       minobs <- -1;  
       mindis <- MAXVAL;
       for (j in 1:m) {
           if (flag[j] == 1) {
              if (mmsmmsdiss$mmdiss[j] < mindis) {
                  mindis <- mmsmmsdiss$mmdiss[j]
                  minobs <- j
              }
           }
       }
       # find agglomerands clus1 and clus2, with former < latter
       if (minobs < mmsmmsdiss$mm[minobs]) {
          clus1 <- minobs
          clus2 <- mmsmmsdiss$mm[minobs]
       }
       if (minobs > mmsmmsdiss$mm[minobs]) {
          clus2 <- minobs
          clus1 <- mmsmmsdiss$nn[minobs]
       }
       # So, agglomeration of pair clus1 < clus2 defines cluster ncl

       #------------------------------------ Block for subnode labels 
       a[ncl] <- clus1                       # aine, or left child node
       b[ncl] <- clus2                       # benjamin, or right child node
       # Now build up ia, ib as version of a, b which is R-compliant
       if (card[clus1] == 1) ia[ncl] <- (-clus1)     # singleton
       if (card[clus2] == 1) ib[ncl] <- (-clus2)     # singleton
       if (card[clus1] > 1) {                # left child is non-singleton
          lastind <- 0
          for (j2 in (m-1):(ncl+1)) {        # Must have m-1 >= ncl+1 here
             if (a[j2] == clus1) lastind <- j2    # Only concerns a[j2]
          }
          ia[ncl] <- m - lastind             # label of non-singleton
       }
       if (card[clus2] > 1) {                # right child is non-singleton
          lastind <- 0
          for (j2 in (m-1):(ncl+1)) {        # Must have m-1 >= ncl+1 here
             if (a[j2] == clus2) lastind <- j2    # Can only concern a[j2]
          }
          ib[ncl] <- m - lastind             # label of non-singleton
       }
       if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
          left <- min(ia[ncl],ib[ncl])
          right <- max(ia[ncl],ib[ncl])
          ia[ncl] <- left                    # Just get left < right
          ib[ncl] <- right
       }
       #--------------------------------------------------------------------

       lev[ncl] <- mindis
       for (j in 1:m) {
           clusmat[j,ncl] <- clusmat[j,ncl+1]
           if (clusmat[j,ncl] == clus2) clusmat[j,ncl] <- clus1
       }

       # store variance decomposition for AMCLUST  PROBABLY NOT NECESSARY... IN lev
#       var.decomp[steps] <- mindis

       # Next we need to update diss array
#       steps <- steps + 1
       # First, the amalgamation
       data.temp <- data.pro2
       for(j in 1:m) {
         foo <- which(clusmat[,ncl]==j)
         if(length(foo)>1) {
           data.temp[,foo[1]] <- apply(data.pro2[,foo], 1, sum, na.rm=TRUE)
           data.temp[,foo[-1]] <- NA
           diss[foo[-1],] <- 99
           diss[,foo[-1]] <- 99
         }
       }
         data.pro2 <- data.temp
  #       mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
         card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
         # Cluster label clus2 is knocked out; following not nec. but no harm
         flag[clus2] <- 0
         mmsmmsdiss$mmdiss[clus2] <- MAXVAL
         for (j1 in 2:m) {
           if ( (j1 != clus2) && (flag[j1] == 1) ) {
             for (j2 in 1:(j1-1)) {
               if ( (j2 != clus2) && (flag[j2] == 1) ) {
  # amalgamate the parts in data.pro (which are in columns)
                 data.foo      <- data.pro2
                 foo.j1j2      <- as.numeric(data.pro2[,j1] + data.pro2[,j2])
                 data.foo[,j1] <- foo.j1j2
                 exclude       <- c(j2, which(flag==0))
                 if(ncl>2) {
#                   if(length(weight==1) & !weight[1]) foo.rda  <- rda(data.clr ~ CLR(data.foo[,-exclude], weight=FALSE)$LR)
#                   if(length(weight==1) & weight[1])  foo.rda  <- rda(data.clr ~ CLR(data.foo[,-exclude], weight=TRUE)$LR)
#                   if(length(weight) > 1) foo.rda  <- rda(data.clr ~ CLR(data.foo[,-exclude], weight=FALSE)$LR)
                   if(length(weight==1) & !weight[1]) foo.RDA  <- RDA(data.clr, cov=CLR(data.foo[,-exclude], weight=FALSE)$LR, weight=cm)
                   if(length(weight==1) & weight[1])  foo.RDA  <- RDA(data.clr, cov=CLR(data.foo[,-exclude], weight=TRUE)$LR, weight=cm)
                   if(length(weight) > 1)             foo.RDA  <- RDA(data.clr, cov=CLR(data.foo[,-exclude], weight=TRUE)$LR, weight=cm)
                   diss[j1,j2]   <- diss[j2,j1] <- TOTVAR - sum(foo.RDA$sv^2)
                 }
#                if(ncl>2) foo.rda <- RDA(data.clr, cov=CLR(data.foo[,c(1,5,9)], weight=FALSE)$LR, weight=FALSE)
#                if(ncl<=2) foo.rda <- rda(data.clr ~ 1)
# NOTE THAT $CA is the unexplained part, which is what you want to have in lev
#                diss[j1,j2]   <- diss[j2,j1] <- (foo.rda$CA$tot.chi/m)*(n-1)/n

                 if(ncl<=2) {
#                   if(length(weight==1) & !weight[1]) diss[j1,j2] <- diss[j2,j1] <- LR.VAR(CLR(CLOSE(data.pro), weight=FALSE))
#                   if(length(weight)==1 & weight[1])  diss[j1,j2] <- diss[j2,j1] <- LR.VAR(CLR(CLOSE(data.pro), weight=TRUE))
#                   if(length(weight) > 1)             diss[j1,j2] <- diss[j2,j1] <- LR.VAR(CLR(CLOSE(data.pro), weight=cm))
                   diss[j1,j2] <- diss[j2,j1] <- LR.VAR(CLR(CLOSE(data.pro), weight=cm))
                 }
               }
             }
          }
        }

#       mass[clus2] <- 0.0
#       for (i in 1:n) {
#           diss[clus2,i] <- MAXVAL
#           diss[i,clus2] <- diss[clus2,i]

       # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
       # i.e. nearest neighbors and the nearest neigh. dissimilarity
#       nnsnnsdiss <- getnns(diss, flag)

### replaced with code to find smallest distance, no error checking needed
### flag all 1s, no need to check

   mm <- rep(0, nrow(diss))
   mmdiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   for (j1 in 1:nrow(diss)) {
          minobs <- -1
          mindis <- MAXVAL
          for (j2 in 1:ncol(diss)) {
              if ( (diss[j1,j2] < mindis) && (j1 != j2) ) {
                 mindis <- diss[j1,j2]
                 minobs <- j2
              }
          }
          mm[j1] <- minobs
          mmdiss[j1] <- mindis
   }

### compatibility with hierclust after "function" call

   mmsmmsdiss <- list(mm=mm, mmdiss=mmdiss)
   }

   temp <- cbind(a,b)
   merge2 <- temp[nrow(temp):1, ]
   temp <- cbind(ia,ib)
   merge <- temp[nrow(temp):1,]
   dimnames(merge) <- NULL
   # merge is R-compliant; later suppress merge2

   #-------------------------------- Build R-compatible order from ia, ib
   orderlist <- c(merge[m-1,1], merge[m-1,2])
   norderlist <- 2
   for (i in 1:(m-2)) {           # For precisely n-2 further node expansions
       for (i2 in 1:norderlist) {       # Scan orderlist
           if (orderlist[i2] > 0) {     # Non-singleton to be expanded
              tobeexp <- orderlist[i2]
              if (i2 == 1) {
                 orderlist <- c(merge[tobeexp,1],merge[tobeexp,2],
                                orderlist[2:norderlist])
              }
              if (i2 == norderlist) {
                 orderlist <- c(orderlist[1:(norderlist-1)],
                                merge[tobeexp,1],merge[tobeexp,2])
              }
              if (i2 > 1 && i2 < norderlist) {
                 orderlist <- c(orderlist[1:(i2-1)], 
                                merge[tobeexp,1],merge[tobeexp,2],
                                orderlist[(i2+1):norderlist])
              }
              norderlist <- length(orderlist)
           }
        }
   }
   orderlist <- (-orderlist)
   class(orderlist) <- "integer              "
   
   xcall <- "hierclust(a,wt)"
   class(xcall) <- "call"

#   retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
#         labels=rownames(a),method="ward",
#         dist.method="euclidean")
#  output square roots of heights, so that sum of squares gives TotVar
   retlist <- list(merge=merge,height=lev[(m-1):1],order=orderlist,
                   method="amalgamation", labels=colnames(data))
   class(retlist) <- "hclust"
   retlist
}