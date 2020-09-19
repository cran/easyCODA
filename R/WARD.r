WARD <- function(LRdata, weight = TRUE, row.wt = NA) {

### easyCODA preparation 

# should be a matrix of logratios or a list object with $LR
# check introduced 16/12/2018
    
    if (is.data.frame(LRdata) | is.matrix(LRdata)) {
        LRfoo <- as.matrix(LRdata)
        if(abs(sum(LRfoo)-nrow(LRfoo))<0.001 | abs(sum(LRfoo)- 100*nrow(LRfoo))<0.1 ) 
            stop("This seems to be compositional matrix, should be logratios")
        weights <- rep(1/ncol(LRfoo), ncol(LRfoo))
    }
    if (is.list(LRdata) & !is.data.frame(LRdata)) {
       LRfoo <- LRdata$LR
       if(abs(sum(LRfoo)-nrow(LRfoo))<0.001 | abs(sum(LRfoo)- 100*nrow(LRfoo))<0.1 ) 
            stop("This seems to be compositional matrix, should be logratios")

        if (!weight[1]) 
            weights <- rep(1/ncol(LRfoo), ncol(LRfoo))
        if (weight[1] & (ncol(LRfoo) > 1)) 
            weights <- LRdata$LR.wt
        if (length(weight) == ncol(LRfoo)) {
            if (sum(weight <= 0) > 0) 
                stop("Some weights zero or negative")
            if (sum(weight) != 1) 
                print("Sum of column weights not exactly 1, but are rescaled")
            weights <- weight/sum(weight)
        }
    }
    if (is.na(row.wt[1])) 
        row.wt <- rep(1/nrow(LRfoo), nrow(LRfoo))
    if (length(row.wt) != nrow(LRfoo)) 
        stop("Error: row weights not the same number as rows of data")
    if (sum(row.wt) != 1) 
        print("Sum of row weights not exactly 1, but are rescaled")
    row.wt <- row.wt/sum(row.wt)

### make equivalent to hierclust

    a <- LRfoo
    wt <- row.wt


   n <- nrow(a)    
   m <- ncol(a)                           
#   diss <- dissim(a, wt)                      # call to function dissim

### replaced by easyCODA computation of weighted logratio distances, squared!

   diss <- matrix(0, n, n)
   rownames(diss) <- colnames(diss) <- rownames(a)

   for (i1 in 2:n) {
       for (i2 in 1:(i1-1)) {
           temp <- 0.0
           for (j in 1:m) {
               # We use the squared Euclidean distance, weighted row & column wise
               temp <- temp + (wt[i1]*wt[i2])/(wt[i1]+wt[i2]) *
                    weights[j] * (a[i1,j]-a[i2,j])^2 
           }
           diss[i2,i1] <- diss[i1,i2] <- temp
       }
   }

   flag <- rep(1, n)                          # active/dead indicator
   a <- rep(0, n-1)                           # left subnode on clustering
   b <- rep(0, n-1)                           # right subnode on clustering
   ia <- rep(0, n-1)                          # R-compatible version of a
   ib <- rep(0, n-1)                          # R-compatible version of b
   lev <- rep(0, n-1)                         # level or criterion values
   card <- rep(1, n)                          # cardinalities
   mass <- wt
   order <- rep(0, n)                         # R-compatible order for plotting

#   nnsnnsdiss <- getnns(diss, flag)           # call to function getnns

### replaced with code to find smallest distance, no error checking needed
### flag all 1s, no need to check

   nn <- rep(0, nrow(diss))
   nndiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   for (i1 in 1:nrow(diss)) {
          minobs <- -1
          mindis <- MAXVAL
          for (i2 in 1:ncol(diss)) {
              if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
                 mindis <- diss[i1,i2]
                 minobs <- i2
              }
          }
          nn[i1] <- minobs
          nndiss[i1] <- mindis
   }

### compatibility with hierclust after "function" call

   nnsnnsdiss <- list(nn=nn, nndiss=nndiss)



   clusmat <- matrix(0, n, n)                 # cluster memberships
   for (i in 1:n) clusmat[i,n] <- i           # init. trivial partition

   MAXVAL <- 1.0e12

   for (ncl in (n-1):1) {                      # main loop 
       # check for agglomerable pair
       minobs <- -1;  
       mindis <- MAXVAL;
       for (i in 1:n) {
           if (flag[i] == 1) {
              if (nnsnnsdiss$nndiss[i] < mindis) {
                  mindis <- nnsnnsdiss$nndiss[i]
                  minobs <- i
              }
           }
       }
       # find agglomerands clus1 and clus2, with former < latter
       if (minobs < nnsnnsdiss$nn[minobs]) {
          clus1 <- minobs
          clus2 <- nnsnnsdiss$nn[minobs]
       }
       if (minobs > nnsnnsdiss$nn[minobs]) {
          clus2 <- minobs
          clus1 <- nnsnnsdiss$nn[minobs]
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
          for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
              if (a[i2] == clus1) lastind <- i2    # Only concerns a[i2]
          }
          ia[ncl] <- n - lastind             # label of non-singleton
       }
       if (card[clus2] > 1) {                # right child is non-singleton
          lastind <- 0
          for (i2 in (n-1):(ncl+1)) {        # Must have n-1 >= ncl+1 here
              if (a[i2] == clus2) lastind <- i2    # Can only concern a[i2]
          }
          ib[ncl] <- n - lastind             # label of non-singleton
       }
       if (ia[ncl] > 0 || ib[ncl] > 0) {     # Check that left < right
          left <- min(ia[ncl],ib[ncl])
          right <- max(ia[ncl],ib[ncl])
          ia[ncl] <- left                    # Just get left < right
          ib[ncl] <- right
       }
       #--------------------------------------------------------------------

       lev[ncl] <- mindis
       for (i in 1:n) {
           clusmat[i,ncl] <- clusmat[i,ncl+1]
           if (clusmat[i,ncl] == clus2) clusmat[i,ncl] <- clus1
       }
       # Next we need to update diss array
       for (i in 1:n) {
           if ( (i != clus1) && (i != clus2) && (flag[i] == 1) ) {
              diss[clus1,i] <- 
      ((mass[clus1]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,i] +
      ((mass[clus2]+mass[i])/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus2,i] -
      (mass[i]/(mass[clus1]+mass[clus2]+mass[i]))*diss[clus1,clus2] 
              diss[i,clus1] <- diss[clus1,i]
           }
       }
       mass[clus1] <- mass[clus1] + mass[clus2]    # Update mass of new cluster
       card[clus1] <- card[clus1] + card[clus2]    # Update card of new cluster
       # Cluster label clus2 is knocked out; following not nec. but no harm
       flag[clus2] <- 0
       nnsnnsdiss$nndiss[clus2] <- MAXVAL
       mass[clus2] <- 0.0
       for (i in 1:n) {
           diss[clus2,i] <- MAXVAL
           diss[i,clus2] <- diss[clus2,i]
       }
       # Finally update nnsnnsdiss$nn and nnsnnsdiss$nndiss
       # i.e. nearest neighbors and the nearest neigh. dissimilarity
#       nnsnnsdiss <- getnns(diss, flag)

### replaced with code to find smallest distance, no error checking needed
### flag all 1s, no need to check

   nn <- rep(0, nrow(diss))
   nndiss <- rep(0.0, nrow(diss))
   MAXVAL <- 1.0e12
   for (i1 in 1:nrow(diss)) {
          minobs <- -1
          mindis <- MAXVAL
          for (i2 in 1:ncol(diss)) {
              if ( (diss[i1,i2] < mindis) && (i1 != i2) ) {
                 mindis <- diss[i1,i2]
                 minobs <- i2
              }
          }
          nn[i1] <- minobs
          nndiss[i1] <- mindis
   }

### compatibility with hierclust after "function" call

   nnsnnsdiss <- list(nn=nn, nndiss=nndiss)

   }

   temp <- cbind(a,b)
   merge2 <- temp[nrow(temp):1, ]
   temp <- cbind(ia,ib)
   merge <- temp[nrow(temp):1,]
   dimnames(merge) <- NULL
   # merge is R-compliant; later suppress merge2

   #-------------------------------- Build R-compatible order from ia, ib
   orderlist <- c(merge[n-1,1], merge[n-1,2])
   norderlist <- 2
   for (i in 1:(n-2)) {           # For precisely n-2 further node expansions
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
   class(orderlist) <- "integer"
   
   xcall <- "hierclust(a,wt)"
   class(xcall) <- "call"

#   retlist <- list(merge=merge,height=as.single(lev[(n-1):1]),order=orderlist,
#         labels=rownames(a),method="ward",
#         dist.method="euclidean")
#  output square roots of heights, so that sum of squares gives TotVar
   retlist <- list(merge=merge,height=sqrt(lev[(n-1):1]),order=orderlist)
   class(retlist) <- "hclust"
   retlist
}