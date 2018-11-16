PLOT.PCA <- function(obj, map="symmetric", rescale=1, dim=c(1,2), main="", 
                     cex=c(0.8,0.8), col=c("blue","red"), colarrows="pink",
                     axes.inv=c(1,1)) {

# updated to fix up colour definitions and character sizes for rows and columns

# plotting function for PCA objects
# map = "symmetric" / "asymmetric" / "contribution"
# rescale = rescaling of column coordinates
# dim = dimensions for plotting, default c(1,2)
# axes.inv = c(1,1) (default); e.g. = c(1,-1) for inverting second axis, or c(-1,-1) for both

  obj.rpc <- obj$rowcoord[,dim] %*% diag(obj$sv[dim] * axes.inv)
  obj.csc <- obj$colcoord[,dim] %*% diag(axes.inv)
  obj.cpc <- obj.csc %*% diag(obj$sv[dim])
  obj.ccc <- obj.csc * sqrt(obj$colmass)

  if(map == "symmetric")    obj.crd <- obj.cpc
  if(map == "asymmetric")   obj.crd <- obj.csc
  if(map == "contribution") obj.crd <- obj.ccc

  perc.hor <- 100 * obj$sv[dim[1]]^2 / sum(obj$sv^2)
  perc.ver <- 100 * obj$sv[dim[2]]^2 / sum(obj$sv^2)

  if(rescale == 1) {
    plot(1.05*rbind(obj.rpc, rescale*obj.crd), type="n", asp=1, 
         xlab=paste("PCA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("PCA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         main=main)
  }
  if(rescale !=1) {
    plot(1.05*rbind(obj.rpc, rescale*obj.crd), type="n", asp=1, 
         xlab=paste("PCA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("PCA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         xaxt="n", yaxt="n", main=main)
    axis(1)
    axis(2)
    axis(3, at=axTicks(3), labels=round(axTicks(3)/rescale, 1), col="black", col.ticks=col[2], col.axis=col[2])
    axis(4, at=axTicks(4), labels=round(axTicks(4)/rescale, 1), col="black", col.ticks=col[2], col.axis=col[2])
  }
  abline(h=0, v=0, col="gray", lty=2)
  if(map != "symmetric") arrows(0, 0, 0.95*rescale*obj.crd[,1], 0.95*rescale*obj.crd[,2], length=0.1, angle=10, col="pink", lwd=2)
  text(obj.rpc, labels=obj$rownames, col=col[1], font=2, cex=cex[1])
  text(rescale*obj.crd, labels=obj$colnames, col=col[2], cex=cex[2], font=4)
}

