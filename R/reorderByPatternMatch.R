#' Reorder By Pattern Match
#'
#' @param P matrix to be matched
#' @param matchTo matrix to match P to
#' @return matched patterns
reorderByPatternMatch <- function(P, matchTo)
{
    # check that P and the matchTo matrix have the same dimensions
    if (nrow(matchTo) != nrow(P) | ncol(matchTo) != ncol(P))
    {
        stop('dimensions of P and matchTo must agree')
    }

    # ensuring that rownames match for simplicty of matching process
    row.names(matchTo) <- row.names(P)

    # compute the correlation between each entry
    corP <- cor(t(matchTo),t(P))

    # initalize the new matrix
    pMatch <- rep(0, nrow(P))
    names(pMatch) <- row.names(P)

    # match patterns in order of correlation
    for (p in 1:(nrow(P)-1))
    {
        ptemp <- which(corP==max(corP),arr.ind=TRUE)
        pMatch[row.names(corP)[ptemp[1]]] <- colnames(corP)[ptemp[2]]
        if (length(corP) > 1)
        {
            corP <- corP[-ptemp[1],-ptemp[2]]
        }
    }
    pMatch[which(pMatch==0)] <- setdiff(names(pMatch), pMatch)

    # return matched patterns
    return(P[pMatch,])
}
