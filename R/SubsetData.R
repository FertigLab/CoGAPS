#' use user provided subsets
#' @keywords internal
#'
#' @param allParams list of all CoGAPS parameters
#' @param total total number of rows (cols) that are being paritioned
#' @return list of subsets
sampleWithExplictSets <- function(allParams, total)
{
    if (all(sapply(allParams$gaps@explicitSets, function(s) is(s, "numeric"))))
    {
        gapsCat(allParams, "using provided indexed subsets\n")
        return(allParams$gaps@explicitSets)
    }
    else if (all(sapply(allParams$gaps@explicitSets, function(s) is(s, "character"))))
    {
        gapsCat(allParams, "using provided named subsets\n")
        if (allParams$gaps@distributed == "genome-wide")
            allNames <- allParams$geneNames
        else
            allNames <- allParams$sampleNames
        return(lapply(allParams$gaps@explicitSets, function(set) which(allNames %in% set)))
    }
}

#' subset rows (cols) proportional to the user provided weights
#' @keywords internal
#'
#' @param allParams list of all CoGAPS parameters
#' @param setSize the size of each subset of the total
#' @return list of subsets
sampleWithAnnotationWeights <- function(allParams, setSize)
{
    # sort annotation group and weights so they match up
    weight <- allParams$gaps@samplingWeight
    weight <- weight[order(names(weight))]
    groups <- unique(allParams$gaps@samplingAnnotation)
    groups <- sort(groups)

    # sample accordingly
    return(lapply(1:allParams$gaps@nSets, function(i)
    {
        groupCount <- sample(groups, size=setSize, replace=TRUE, prob=weight)
        sort(unlist(sapply(groups, function(g)
        {
            groupNdx <- which(allParams$gaps@samplingAnnotation == g)
            sample(groupNdx, size=sum(groupCount == g), replace=TRUE)
        })))
    }))
}

#' subset data by uniformly partioning rows (cols)
#' @keywords internal
#'
#' @param allParams list of all CoGAPS parameters
#' @param total total number of rows (cols) that are being paritioned
#' @param setSize the size of each subset of the total
#' @return list of subsets
sampleUniformly <- function(allParams, total, setSize)
{
    sets <- list()
    remaining <- 1:total
    for (n in 1:(allParams$gaps@nSets - 1))
    {
        selected <- sample(remaining, setSize, replace=FALSE)
        sets[[n]] <- sort(selected)
        remaining <- setdiff(remaining, selected)
    }
    sets[[allParams$gaps@nSets]] <- sort(remaining)
    return(sets)
}

#' partition genes/samples into subsets
#' @keywords internal
#'
#' @description either genes or samples or partitioned depending on the type
#' of distributed CoGAPS (i.e. genome-wide or single-cell)
#' @param data either file name or matrix
#' @param allParams list of all CoGAPS parameters
#' @return list of sorted subsets of either genes or samples
createSets <- function(data, allParams)
{
    subsetRows <- xor(allParams$transposeData, allParams$gaps@distributed == "genome-wide")
    total <- ifelse(subsetRows, nrowHelper(data), ncolHelper(data))
    setSize <- floor(total / allParams$gaps@nSets)

    gapsCat(allParams, "Creating subsets...")

    if (!is.null(allParams$gaps@explicitSets))
    {
        if (length(allParams$gaps@explicitSets) != allParams$gaps@nSets)
            stop("nSets does not match number of explicit sets given")
        return(sampleWithExplictSets(allParams, total))
    }

    if (!is.null(allParams$gaps@samplingAnnotation))
    {
        gapsCat(allParams, "sampling with annotation weights\n")
        return(sampleWithAnnotationWeights(allParams, setSize))
    }
    gapsCat(allParams, "\n")
    return(sampleUniformly(allParams, total, setSize))
}
