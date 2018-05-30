invALR <- function(ALRmatrix, part.names=paste("part",1:(ncol(ALRmatrix)+1),sep=""), denom=NA) {
# Inverts a I x (J-1) matrix of additive logratios (ALRs)
# back to their original compositional part values
# ALRmatrix = matrix of ALRs, ratios in columns
# Returns an I x J matrix of parts, where last column is the reference part 
  parts <- exp(ALRmatrix) / (1 + apply(exp(ALRmatrix), 1, sum))
  parts <- cbind(parts, 1-apply(parts, 1, sum))
  colnames(parts) <- part.names
# rearrange if denominator specified
  if(!is.na(denom)) {
    foo <- parts
    foo[,denom] <- parts[,ncol(parts)]
    foo[,(1:ncol(parts))[-denom]] <- parts[, 1:(ncol(parts)-1)]
    parts <- foo
  }
  parts
}
