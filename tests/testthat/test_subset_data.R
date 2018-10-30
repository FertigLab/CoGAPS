context("CoGAPS")

getIndex <- function(s)
    as.numeric(strsplit(s, "_")[[1]][2])

getIndices <- function(set)
    unname(sapply(set, getIndex))

countType <- function(set, type, anno)
    sum(getIndices(set) %in% which(anno == type))

getHistogram <- function(set, anno)
    sapply(letters[1:5], function(letter) countType(set, letter, anno))

test_that("Subsetting Data",
{
    data(GIST)
    data(SimpSim)
    testMatrix <- GIST.matrix

    # sampling with weights
    anno <- sample(letters[1:5], size=nrow(testMatrix), replace=TRUE)
    weights <- c(0,1,4,1,0)
    names(weights) <- letters[1:5]
    params <- new("CogapsParams")
    params <- setAnnotationWeights(params, annotation=anno, weights=weights)
    result <- GWCoGAPS(testMatrix, params, messages=FALSE, seed=123,
        nIterations=5000, geneNames=paste("Gene", 1:nrow(testMatrix), sep="_"))

    hist <- sapply(getSubsets(result), getHistogram, anno=anno)
    freq <- unname(rowSums(hist) / sum(hist))
    
    expect_true(all.equal(freq, unname(weights / sum(weights)), tol=0.1))

    # running cogaps with given subsets
    sets <- list(1:225, 226:450, 451:675, 676:900)
    result <- GWCoGAPS(SimpSim.data, nPatterns=7, explicitSets=sets,
        messages=FALSE, matchedPatterns=matrix(1,nrow=ncol(SimpSim.data),ncol=7))
    subsets <- lapply(getSubsets(result), getIndices)
    expect_true(all(sapply(1:4, function(i) all.equal(sets[[i]], subsets[[i]]))))
})