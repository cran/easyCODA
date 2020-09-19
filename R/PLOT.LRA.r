PLOT.LRA <- function(obj, map="symmetric", rescale=1, dim=c(1,2), axes.inv = c(1,1), main="", cols=c("blue","red"), colarrows = "pink", cexs=c(0.8,0.8), fonts=c(2,4)) {
# plotting function for LRA objects
# obj       = LRA object
# map       = "symmetric" / "asymmetric" / "contribution"
# rescale   = rescaling of column coordinates
# dim       = selected dimensions
# axes.inv  = c(1,1) (default); e.g. = c(1,-1) for inverting second axis, or c(-1,-1) for both
# main      =  main title
# cols      = c("blue","red") (default) colours for row and column labels
# colarrows = "pink" (default) colour for arrows in asymmetric and contribution biplots
# cexs      = c(0.8,0.8) (default) character explansion factors for row and column labels
# fonts     = c(2,4) (default) font choice for rows and column labels)

  obj.rpc <- obj$rowcoord[,dim] %*% diag(obj$sv[dim] * axes.inv)
  obj.csc <- obj$colcoord[,dim] %*% diag(axes.inv)
  obj.cpc <- obj.csc %*% diag(obj$sv[dim])
  obj.ccc <- obj.csc * sqrt(obj$colmass)

  if(map == "symmetric")    obj.crd <- obj.cpc
  if(map == "asymmetric")   obj.crd <- obj.csc
  if(map == "contribution") obj.crd <- obj.ccc
  if((map != "symmetric") & (map != "asymmetric") & (map != "contribution"))
     stop("map option is not a valid choice of present options: symmetric, asymmetric, contribution")

  perc.hor <- 100 * obj$sv[dim[1]]^2 / sum(obj$sv^2)
  perc.ver <- 100 * obj$sv[dim[2]]^2 / sum(obj$sv^2)

  if(rescale == 1) {
    plot(1.05*rbind(obj.rpc, rescale*obj.crd), type="n", asp=1, 
         xlab=paste("LRA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("LRA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         main=main)
  }
  if(rescale !=1) {
    plot(1.05*rbind(obj.rpc, rescale*obj.crd), type="n", asp=1, 
         xlab=paste("LRA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("LRA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         xaxt="n", yaxt="n", main=main)
    axis(1)
    axis(2)
    axis(3, at=axTicks(3), labels=round(axTicks(3)/rescale,2), col=cols[2], col.ticks=cols[2], col.axis=cols[2])
    axis(4, at=axTicks(4), labels=round(axTicks(4)/rescale,2), col=cols[2], col.ticks=cols[2], col.axis=cols[2])
  }
  abline(h=0, v=0, col="gray", lty=2)
  if(map != "symmetric") arrows(0, 0, 0.95*rescale*obj.crd[,1], 0.95*rescale*obj.crd[,2], length=0.1, angle=10, col=colarrows, lwd=2)
  text(obj.rpc, labels=obj$rownames, col=cols[1], font=fonts[1], cex=cexs[1])
  text(rescale*obj.crd, labels=obj$colnames, col=cols[2], font=fonts[2], cex=cexs[2])
}


