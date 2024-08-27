CA <- function(data, nd=2, suprow = NA, supcol = NA) {

# This is mainly a wrapper for the ca function with some alternative options to come

# CA computations
# just make ca object 
  data.ca             <- ca(data, suprow=suprow, supcol=supcol)  
# above caold is temporary fix until vegan updated
  data.ca$rowpcoord   <- data.ca$rowcoord[,1:nd] %*% diag(data.ca$sv[1:nd])
  data.ca$colpcoord   <- data.ca$colcoord[,1:nd] %*% diag(data.ca$sv[1:nd])
# keep class as "ca", not "CA"
  class(data.ca)      <- "ca"
  return(data.ca)
  }

