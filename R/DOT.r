DOT <- function(data, cols=NA, names=NA, groups=NA, pch=NA, horizon=FALSE, jitter=1, 
                xscale=NA, xscalefac=1, yaxis=TRUE, shownames=TRUE, main="", ylab="", 
                xlim=c(NA,NA), ylim=c(NA, NA), cex=1)
{

# function to plot dots, assuming diffe	rent samples in columns of data, or defined by 
# group code.  A common scale is assumed.

# data: matrix with data in columns (even if single sample, should be a matrix)
# cols: colours of points for each sample, default rainbow
# names: labels for variables, by default the column names of data, or group names
# groups: group codes to split the data vector into separate plots
# main: title for figure
# ylab: y-axis label
# ylim: user-defined y limits, otherwise computed by function
# yaxis: =F if no y-axis plotted, by default=T
# cex: expansion factor for symbols
# pch: alternative set of symbols for samples, by default small circles (pch=1)
# horizon: T if horizontal gray dashed lines required at "nice" y-values, not implemented yet
# jitter: 1 by default, increase or decrease slightly for more jitter
# xscalefac: =1, can be changed to change gaps between plots
# shownames: =T, =F to not show group names and add them externally

# determine x- and y-ranges
  data <- as.matrix(data)
  if(is.na(xlim[1])) {
    xrange <- c(0, ncol(data)+1)
    xscale <- xscalefac*(1:ncol(data))
  }
  if(!is.na(groups[1])) {
    ngroups <- length(unique(groups))
    xrange <- c(0, ngroups+1)
    xscale <- xscalefac*(1:ngroups)
  }
  if(!is.na(xlim[1]) & !is.na(xlim[2]) & !is.na(xscale[1])) {
    xrange <- xlim
    xscale <- xscale
  }
  if(is.na(ylim[1])) {
    yrange <- max(data, na.rm=T)-min(data, na.rm=T)
    ylim   <- c(min(data, na.rm=T), max(data, na.rm=T)+0.05*yrange)
  }

# DATA is the matrix of values, by default it is = data

  if(is.na(groups[1])) {
    DATA <- data
    if(is.na(names[1])) names <- colnames(data)
  }
  if(!is.na(groups[1])) {
    DATA <- matrix(NA, nrow=max(table(groups)), ncol=ngroups)
    colnames(DATA) <- unique(groups)
    for(j in 1:ncol(DATA)) {
      foo <- data[groups==colnames(DATA)[j]] 
      DATA[1:length(foo),j] <- foo
    }
    names <- colnames(DATA)
    if(is.na(cols[1])) cols <- rainbow(ncol(DATA))
    if(length(pch)==1) pch <- rep(19,ncol(DATA))
  }
# plot
  if(is.na(cols[1])) cols <- rainbow(ncol(DATA))
  if(is.na(pch[1])) pch <- rep(19,ncol(DATA))
  par(mar=c(4,4.2,3,0), font.lab=2, las=1, cex.lab=1.2)
  plot(0, 0, type="n", xlim=xrange, ylim=ylim, xaxt="n", yaxt="n", bty="n", 
       xlab="", main=main, ylab=ylab)
  if(yaxis) axis(2)
  axis(1, at=xscale, labels=names, tick=F, font=2, col=cols)
  
  for(j in 1:ncol(DATA)) {
    foodata <- DATA[!is.na(DATA[,j]),j]
    footab  <- table(foodata)
    yval    <- as.numeric(names(footab))
    for(i in 1:length(yval)) {
      if(footab[[i]]==1) points(xscale[j],yval[i], pch=pch[j], col=cols[j], cex=cex)
      if(footab[[i]]>1) {
        foo <- footab[[i]]
        foo <- c(1:foo)-mean(1:foo)
        foo <- foo * 0.05 * jitter
        foo <- xscale[j] + foo
        points(foo, rep(yval[i], length(foo)), pch=pch[j], col=cols[j], cex=cex, font=2)
      } 
    }
  }
}

