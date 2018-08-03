PLOT.RDA <- function(obj, map="symmetric", indcat=NA, rescale=1, dim=c(1,2), main="", 
                     axes.inv=c(1,1), rowstyle=1, cols=c("blue","red","forestgreen"), 
                     colarrows=c("pink","lightgreen"), colrows=NA, pchrows=NA, colcat=NA,
                     cexs=c(0.8,0.8,0.8), fonts=c(2,4,4)) {
# plotting function for RDA objects
# obj        RDA object 
# indcat     positions of dummy or fuzzy variables in covariates
# map        "symmetric" / "asymmetric" / "contribution"
# rescale    rescaling of column coordinates and covariate regression coefficients
# axes.inv   inversion of axes (default: none); e.g. = c(1,-1) for inverting second axis, 
#            or c(-1,-1) for both
# rowstyle   = 1 (default), = 2 (samples as supplementary points)
# cols       colours for row, column and covariate labels (default: c("blue","red","forestgreen"))
# colarrows  colour for arrows in asymmetric and contribution biplots
# colrows    (optional) vector of colours for rows
# pchrows    (optional) vector of symbols for rows
# colcat     (optional) vector of colours for category means
# cexs       expansion/contraction factors for row, column and covariate labels (default: c(0.8,0.8,0.8))
# fonts      fonts for row, column and covariate labels (default: c(2,4,4))

  if(length(rescale)==1) rescale <- c(rescale, rescale)
  vars <- obj$cov
  var.ind <- 1:ncol(vars)
  if(ncol(vars)==length(indcat)) var.ind <- NA
  if((ncol(vars)>length(indcat)) & !is.na(indcat[1])) var.ind <- (1:ncol(vars))[-indcat]

  if(rowstyle==1) obj.rpc <- obj$rowpcoord[,dim] %*% diag(axes.inv)
  obj.csc <- obj$colcoord[,dim] %*% diag(axes.inv)

  if(rowstyle==2) {
    data.centred <- scale(obj$N, scale=FALSE)
    obj.rpc <- data.centred %*% diag(sqrt(obj$colmass)) %*% obj.csc
  }

  obj.cpc <- obj$colpcoord[,dim] %*% diag(axes.inv)
  obj.ccc <- obj.csc * sqrt(obj$colmass)

  obj.cvcrd <- obj$covcoord[,dim] %*% diag(axes.inv)

  if(map == "symmetric")    obj.crd <- obj.cpc
  if(map == "asymmetric")   obj.crd <- obj.csc
  if(map == "contribution") obj.crd <- obj.ccc

  perc.hor <- 100 * obj$sv[dim[1]]^2 / sum(obj$sv^2)
  perc.ver <- 100 * obj$sv[dim[2]]^2 / sum(obj$sv^2)

  if(rescale[1] == 1) {
    plot(1.05*rbind(obj.rpc, rescale[1]*obj.crd), type="n", asp=1, 
         xlab=paste("RDA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("RDA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         main=main)
  }
  if(rescale[1] !=1) {
    plot(1.05*rbind(obj.rpc, rescale[1]*obj.crd), type="n", asp=1, 
         xlab=paste("RDA dimension ",dim[1]," (",round(perc.hor,1),"%)", sep=""), 
         ylab=paste("RDA dimension ",dim[2]," (",round(perc.ver,1),"%)", sep=""),
         xaxt="n", yaxt="n", main=main)
    axis(1)
    axis(2)
    axis(3, at=axTicks(3), labels=round(axTicks(3)/rescale[1], 2), col="black", col.ticks="red", col.axis=cols[2])
    axis(4, at=axTicks(4), labels=round(axTicks(4)/rescale[1], 2), col="black", col.ticks="red", col.axis=cols[2])
  }
  if(!is.na(indcat[1])) {
    cov.cat <- (1:nrow(obj$covcoord))[indcat]
    cov.catcoord <- t(as.matrix(obj$cov[,indcat])) %*% obj.rpc / apply(obj$cov[, indcat], 2, sum) 
  }

# rescale continuous covariates
  covmax1 <- max(abs(obj.cvcrd[,1]))
  covmax2 <- max(abs(obj.cvcrd[,2]))
  colmax1 <- max(abs(obj.crd[,1]))
  colmax2 <- max(abs(obj.crd[,2]))
  covscale <- 0.8*max(colmax1/covmax1, colmax2/covmax2)

  abline(h=0, v=0, col="gray", lty=2)
  if(map != "symmetric") arrows(0, 0, 0.95*rescale[1]*obj.crd[,1], 0.95*rescale[1]*obj.crd[,2], length=0.1, angle=10, col=colarrows[1], lwd=2)
  arrows(0, 0, 0.95*rescale[2]*covscale*obj.cvcrd[var.ind,1], 0.95*rescale[2]*covscale*obj.cvcrd[var.ind,2], length=0.1, angle=10, col=colarrows[2], lwd=2)
  if(is.na(colrows[1]) & is.na(pchrows[1])) text(obj.rpc, labels=obj$rownames, col=cols[1], font=fonts[1], cex=cexs[1])
  if(length(colrows)>1 & is.na(pchrows[1])) text(obj.rpc, labels=obj$rownames, col=colrows, font=fonts[1], cex=cexs[1])
  if(is.na(colrows[1]) & length(pchrows)>1) points(obj.rpc, pch=pchrows, col=cols[1], cex=cexs[1])
  if(length(colrows)>1 & length(pchrows)>1) points(obj.rpc, pch=pchrows, col=colrows, cex=cexs[1])
  text(rescale[1]*obj.crd, labels=obj$colnames, col=cols[2], cex=cexs[2], font=fonts[2])
  text(rescale[2]*covscale*obj.cvcrd[var.ind,1],rescale[2]*covscale*obj.cvcrd[var.ind,2], labels=obj$covnames[var.ind], col=cols[3], cex=cexs[3], font=fonts[3])
  if(!is.na(indcat[1]) & is.na(colcat[1])) text(cov.catcoord[,1], cov.catcoord[,2], labels=obj$covnames[indcat], col=cols[3], cex=cexs[3], font=fonts[3])
  if(!is.na(indcat[1]) & !is.na(colcat[1])) text(cov.catcoord[,1], cov.catcoord[,2], labels=obj$covnames[indcat], col=colcat, cex=cexs[3], font=fonts[3])
}

