PLR <- function (data, ordering = 1:ncol(data), weight = TRUE) 
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
    if (sum(intersect(1:ncol(data), ordering) == 1:ncol(data)) != 
        ncol(data)) 
        stop("ordering is not a permutation of the columns of data matrix")
    plr <- matrix(0, nrow(data), ncol(data) - 1)
    data <- data[, ordering]
    weights <- weights[ordering]
    colnames(plr) <- 1:(ncol(data) - 1)
    plr.weights <- rep(0, ncol(data) - 1)
    coeff       <- rep(0, ncol(data) - 1)
    for (j in 1:(ncol(data) - 1)) {
        ilr <- ILR(data, numer = j, denom = (j + 1):ncol(data), 
            weight = weights)
        denom <- rep(0, nrow(data))
        for(k in (j+1):ncol(data)) denom <- denom+weights[k] * log(data[,k])
        coeff[j] <- (weights[j]*sum(weights[(j+1):ncol(data)])) / (weights[j] + sum(weights[(j+1):ncol(data)]))
        plr[, j] <- (log(data[,j]) - denom/sum(weights[(j+1):ncol(data)]))
    }
    plr.weights <- coeff
    colnames(plr) <- names(plr.weights) <- colnames(data)[1:(ncol(data)-1)]
    if (length(rownames(data)) == 0) 
        rownames(plr) <- 1:nrow(data)
    if (length(rownames(data)) > 0) 
        rownames(plr) <- rownames(data)
    return(list(LR = plr, LR.wt = plr.weights))
}
