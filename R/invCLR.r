invCLR <- function(CLRmatrix, part.names=colnames(CLRmatrix)) {
# Inverts a I x J matrix of centred logratios (CLRs)
# back to their original compositional part values
# CLRmatrix = matrix of CLRs, ratios in columns
# Returns an I x J matrix of parts
  clr <- exp(CLRmatrix) / apply(exp(CLRmatrix), 1, sum)
  colnames(clr) <- part.names
  clr
}

