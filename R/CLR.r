CLR <- function (data, weight = TRUE) 
{
    if (sum(data == 0) > 0) 
        stop("Error: some data values are zero")
    data <- as.matrix(data/apply(data, 1, sum))
    if (!weight[1]) 
        weights <- rep(1/ncol(data), ncol(data))
    if (weight[1]) 
        weights <- apply(data, 2, mean)
    if (length(weight) == ncol(data)) {
        if (sum(weight <= 0) > 0) 
            stop("Error: some weights zero or negative")
        if (sum(weight) != 1) 
            print("Sum of weights not exactly 1, but are rescaled")
        weights <- weight/sum(weight)
    }
    clr <- sweep(log(data), 1, apply(log(data) %*% diag(weights), 
        1, sum))
    colnames(clr) <- colnames(data)
    clr.weights <- weights
    names(clr.weights) <- colnames(clr)
    if(length(rownames(data))==0) rownames(clr) <- 1:nrow(data)
    if(length(rownames(data))> 0) rownames(clr) <- rownames(data)   
    return(list(LR = clr, LR.wt = clr.weights))
}
