BAR <- function(data, 
                cols=rainbow(ncol(data)), 
                col.names=colnames(data), 
                row.names=rownames(data),
                order.column=NA, eps=0.5,
                main="", ylab="", 
                ylim=c(0,nrow(data)), xlim=c(0,100), 
                cex=1) {

# order of samples
  neworder <- 1:nrow(data)
  if(!is.na(order.column)) neworder <- order(data[,order.column])
  par(mar=c(4.7,10,3,1), mgp=c(6,0,0), font.lab=2, cex.lab=1.2, xpd=NA)
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
      segments(data.cum[irow.ind,j], irow , data.cum[irow.ind,j+1]-0.5, irow, col=cols[j], lwd=12, lend=1)
    }
  }

# count characters in column names
  col.char <- nchar(col.names)
  tot.char <- 2*sum(col.char) + 3*ncol(data) + 1*ncol(data) + 6*(ncol(data)-1)
  starts <- cumsum(2*col.char+10)[1:2] 
  starts <- c((100-tot.char)/2, (100-tot.char)/2 + starts)
  for(j in 1:ncol(data)) {
    segments(starts[j],-3,starts[j]+3,-3, col=cols[j], lwd=12, lend=1)
  }
  text(starts+4,rep(-3,3), labels=col.names, pos=4, offset=0.2, cex=cex)
}

