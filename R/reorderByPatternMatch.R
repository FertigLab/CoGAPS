# reorderByPatternMatch: function to identify corresponding rows
#                        between two P matrices
# History: EJF - original CoGAPS

# Inputs: P - matrix to be matched
#         matchTo - matrix to match P to

# Output: P with reordered rows to match matchTo

#'\code{reorderByPatternMatch} plots the output A and P matrices as a
#' heatmap and line plot respectively
#'
#'@param P matrix to be matched
#'@param matchTo matrix to match P to
#'@export

reorderByPatternMatch <- function(P, matchTo) {

    # check that P and the matchTo matrix have the same dimensions for valid matching
    if (nrow(matchTo) != nrow(P) | ncol(matchTo) != ncol(P)) {
        stop('CoGAPS: reorderByPatternMatch: dimensions of P and matchTo must agree')
    }

    # ensuring that rownames match for simplicty of matching process
    row.names(matchTo) <- row.names(P)

    # compute the correlation between each entry
    corP <- cor(t(matchTo),t(P))

    # initalize the new matrix
    pMatch <- rep(0, nrow(P))
    names(pMatch) <- row.names(P)

    # match patterns in order of correlation
    for (p in 1:(nrow(P)-1)) {
        ptemp <- which(corP==max(corP),arr.ind=T)
        pMatch[row.names(corP)[ptemp[1]]] <- colnames(corP)[ptemp[2]]
        if (length(corP) > 1) {
            corP <- corP[-ptemp[1],-ptemp[2]]
        }
    }
    pMatch[which(pMatch==0)] <- setdiff(names(pMatch), pMatch)

    # return matched patterns
    return(P[pMatch,])

}
