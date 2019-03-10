LR.VAR <- function (LRdata, row.wt = NA, weight = TRUE, vars = FALSE) 
{   
# should be a matrix of logratios or a list object with $LR
# check introduced 16/12/2018
    
    if (is.data.frame(LRdata) | is.matrix(LRdata)) {
        LRfoo <- as.matrix(LRdata)
        if(abs(sum(LRfoo)-nrow(LRfoo))<0.001 | abs(sum(LRfoo)- 100*nrow(LRfoo))<0.1 ) 
            stop("This seems to be compositional matrix, should be logratios")
        weights <- rep(1/ncol(LRfoo), ncol(LRfoo))
    }
    if (is.list(LRdata) & !is.data.frame(LRdata)) {
       LRfoo <- LRdata$LR
       if(abs(sum(LRfoo)-nrow(LRfoo))<0.001 | abs(sum(LRfoo)- 100*nrow(LRfoo))<0.1 ) 
            stop("This seems to be compositional matrix, should be logratios")

        if (!weight[1]) 
            weights <- rep(1/ncol(LRfoo), ncol(LRfoo))
        if (weight[1] & (ncol(LRfoo) > 1)) 
            weights <- LRdata$LR.wt
        if (length(weight) == ncol(LRfoo)) {
            if (sum(weight <= 0) > 0) 
                stop("Some weights zero or negative")
            if (sum(weight) != 1) 
                print("Sum of column weights not exactly 1, but are rescaled")
            weights <- weight/sum(weight)
        }
    }
    if (is.na(row.wt[1])) 
        row.wt <- rep(1/nrow(LRfoo), nrow(LRfoo))
    if (length(row.wt) != nrow(LRfoo)) 
        stop("Error: row weights not the same number as rows of data")
    if (sum(row.wt) != 1) 
        print("Sum of row weights not exactly 1, but are rescaled")
    row.wt <- row.wt/sum(row.wt)
    if (ncol(LRfoo) == 1) {
        LRmean <- sum(LRfoo * row.wt)
        LRtotvar <- sum((LRfoo - LRmean)^2 * row.wt)
    }
    if (ncol(LRfoo) > 1) {
        LRfoo <- sweep(LRfoo, 2, apply(LRfoo * row.wt, 2, sum))
        LRtotvar <- sum(diag(row.wt) %*% LRfoo^2 %*% diag(weights))
        if (vars) {
            LRvars <- apply(diag(row.wt) %*% LRfoo^2 %*% diag(weights), 
                2, sum)
            names(LRvars) <- colnames(LRfoo)
        }
    }
    if (!vars) 
        return(LRtotvar)
    if (vars) 
        return(list(LRtotvar = LRtotvar, LRvars = LRvars))
}
