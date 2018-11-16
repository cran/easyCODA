CA <- function(data, nd = 2, suprow = NA, supcol = NA) {

# This is mainly a wrapper for the ca function with some alternative options to come

# CA computations
# just make ca object 
  data.ca             <- ca(data, rowsup=suprow, colsup=supcol)  
  data.ca$rowpcoord   <- data.ca$rowcoord %*% diag(data.ca$sv)
  data.ca$colpcoord   <- data.ca$colcoord %*% diag(data.ca$sv)
  return(data.ca)
  }

