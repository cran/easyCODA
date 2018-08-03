DUMMY <- function(x, catnames=NA)

# function to code a factor into dummy variables
# 

{ 
  x.cats <- sort(unique(x))
  if(length(catnames)>1) x.cats <- catnames
  Z <- matrix(0, nrow=length(x), ncol=length(x.cats))
  colnames(Z) <- x.cats
  for(i in 1:length(x)) Z[i,x[i]] <- 1
  Z
}
