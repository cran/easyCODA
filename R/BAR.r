BAR <- function(data, 
                cols=rainbow(ncol(data)), 
                col.names=colnames(data), 
                row.names=rownames(data),
                order.column=NA, eps=0.5,
                main="", ylab="", 
                ylim=c(0,nrow(data)), xlim=c(0,100), 
                cex=1, truncate=NA) {

# warning for large matrices
  if(nrow(data) > 50 | ncol(data) > 15) print("Warning: this plot function only works for smaller matrices (50 rows, 15 columns, at most)")

# first count characters in column names
  if(!is.na(truncate)) col.names <- substr(col.names, 1, truncate)
  col.char <- nchar(col.names)
  spchar <- 20/10
  spbloc <- 12
  tot.char <- spchar*sum(col.char)+spbloc*ncol(data)
  nr.rows  <- ceiling(tot.char/100)
  nr.cols  <- ceiling(length(col.names)/nr.rows)
  starts   <- matrix(0, nr.rows, nr.cols)
  for(row in 1:nr.rows) {
    set.cols <- (nr.cols*(row-1)+1):min((nr.cols*row),length(col.names))
    tot.char.row <- sum(col.char[set.cols]) + (length(set.cols))*spbloc
    border <- (100-tot.char.row)/2
    col1 <- 1
    for(j in 1:length(set.cols)) {
      if(j==1) starts[row,col1] <- border
      if(j>1)  starts[row, col1] <- border + sum(col.char[set.cols[1]:set.cols[j-1]]) + spbloc*(col1-1)
      col1 <- col1+1
    } 
  }
  
# order of samples
  neworder <- 1:nrow(data)
  if(!is.na(order.column)) neworder <- order(data[,order.column])
  longestrowlab <- max(nchar(rownames(data)))
  par(mar=c(2+nrow(data)/15+2*nr.rows,1+longestrowlab/2,2,1), mgp=c(6,0,0), font.lab=2, cex.lab=1.2, xpd=NA)
  plot(0,0,type="n", ylab=ylab, xlab="", main=main, ylim=ylim, xlim=xlim, bty="n", xaxt="n", yaxt="n")
  axis(1, labels=c(0,20,40,60,80,100), at=c(0,20,40,60,80,100), mgp=c(2, 0.7, 0))
  axis(2, tick=FALSE, labels=row.names[neworder], at=seq(nrow(data)-0.5,0.5,-1), las=1, cex.axis=cex)

# loop for plotting, first compute cumulative sums across rows
  data.cum <- cbind(rep(0, nrow(data)),t(apply(data, 1, cumsum)))
  data.cum[, ncol(data.cum)] <- data.cum[, ncol(data.cum)] + eps
  for(i in 1:nrow(data)) {
    irow     <- nrow(data)+0.5-i
    irow.ind <- neworder[i]
    for(j in 1:ncol(data)) {
      segments(data.cum[irow.ind,j], irow , data.cum[irow.ind,j+1]-eps, irow, col=cols[j], lwd=12, lend=1)
    }
  }

# plotting legend
  jlab <- 1
  for(i in 1:nr.rows) {
    for(j in 1:nr.cols) {
      if(starts[i,j]>0) {
        segments(starts[i,j],-1-nrow(data)/15-2.2*i,starts[i,j]+3,-1-nrow(data)/15-2.2*i, col=cols[jlab], lwd=12, lend=1)
        text(starts[i,j]+3,-1-nrow(data)/15-2*i, labels=col.names[jlab], pos=4, offset=0.2, cex=cex)
        jlab <- jlab+1
      }
    }
  }
}

