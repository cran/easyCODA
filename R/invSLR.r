invSLR <- function (SLRmatrix, part.names=NA, ratio.names=colnames(SLRmatrix)) {
# SLRmatrix contains the I x (J-1) matrix of SLRs
# part.names contains the names of the parts (required)
# ratio.names contains the definition of the ratios, using the part names
# SLRdef is the intermediate J-1 x J matrix of 1s and -1s indicating the 
#   numerator and denominator parts

    if(is.na(part.names[1])) stop("Part names not specified")
    parts.given <- NULL
    J <- ncol(SLRmatrix)+1                  # number of parts
    parts <- matrix(0, nrow(SLRmatrix), J)  # resultant inverses
    for(i in 1:nrow(SLRmatrix)) {
      A <- matrix(0, J, J)
      for(j in 1:(J-1)) {
# numerator and denominator of the ratio
        split <- strsplit(ratio.names[j], "/")[[1]]
# numerator parts
        numer <- strsplit(split[1], "&")[[1]]
        if(length(intersect(numer, part.names))!=length(numer)) stop("Unknown part(s)",denom,"in a numerator")
        parts.given <- c(parts.given, numer)
# denominator parts
        denom <- strsplit(split[2], "&")[[1]]
        if(length(intersect(denom, part.names))!=length(denom)) stop(paste("Unknown part(s)",denom,"in a denominator"))
        parts.given <- c(parts.given, denom)
        for(k in 1:length(numer)) {
          A[j, which(part.names==numer[k])]  <- 1
        }
        for(k in 1:length(denom)) {
          A[j, which(part.names==denom[k])] <- -exp(SLRmatrix[i,j])
        }
      }
      if(!all(part.names %in% parts.given)) stop("A part or parts not mentioned in list of SLRs")
      A[J,] <- rep(1, J)
      b <- rep(0, J)
      b[J] <- 1
      parts[i,] <- solve(A, b)
    }
    colnames(parts) <- part.names
    parts
}
