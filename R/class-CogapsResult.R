#' CogapsResult
#' @export
#'
#' @description Contains all output from Cogaps run
#' @importClassesFrom S4Vectors DataFrame Annotated character_OR_NULL
setClass("CogapsResult", contains="Annotated", slots=c(
    sampleFactors = "ANY",      # Pmean transpose
    featureLoadings = "ANY",    # Amean
    featureStdDev = "ANY",      # Asd
    sampleStdDev = "ANY",       # Psd
    NAMES = "character_OR_NULL",
    factorData = "DataFrame"),
)

#' Constructor for CogapsResult
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
function(.Object, Amean, Pmean, Asd, Psd, seed, meanChiSq, diagnostics, ...)
{
    .Object@featureLoadings <- Amean
    .Object@sampleFactors <- t(Pmean)

    if (!is.null(Asd))
        .Object@featureStdDev <- Asd
    if (!is.null(Psd))
        .Object@sampleStdDev <- t(Psd)

    .Object@metadata[["seed"]] <- seed
    .Object@metadata[["meanChiSq"]] <- meanChiSq
    .Object@metadata[["diagnostics"]] <- diagnostics

    .Object <- callNextMethod(.Object, ...)
    .Object
})

setValidity("CogapsResult",
    function(object)
    {
        if (sum(object@featureLoadings < 0) > 0 | sum(object@featureStdDev < 0) > 0)
            "fatal error - negative values in feature x factor Matrix"
        if (sum(object@sampleFactors < 0) > 0 | sum(object@sampleStdDev < 0) > 0)
            "fatal error - negative values in factor x sample Matrix"
    }
)    

################################### GENERICS ###################################

#' return chi-sq of final matrices
#' @export
#' @docType methods
#' @rdname getMeanChiSq-methods
#'
#' @param object an object of type CogapsResult
#' @return chi-sq error
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' meanChiSq(result)
setGeneric("getMeanChiSq", function(object)
    {standardGeneric("getMeanChiSq")})

#' compute z-score matrix
#' @export
#' @docType methods
#' @rdname calcZ-methods
#'
#' @description calculates the Z-score for each element based on input mean
#' and standard deviation matrices
#' @param object an object of type CogapsResult
#' @param which either "feature" or "sample" indicating which matrix to
#' calculate the z-score for
#' @return matrix of z-scores
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' feature_zscore <- calcZ(result)
setGeneric("calcZ", function(object, which="feature")
    {standardGeneric("calcZ")})

#' reconstruct gene
#' @export
#' @docType methods
#' @rdname reconstructGene-methods
#'
#' @param object an object of type CogapsResult
#' @param genes an index of the gene or genes of interest
#' @return the D' estimate of a gene or set of genes
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' D_estimate <- reconstructGene(result)
setGeneric("reconstructGene", function(object, genes=NULL)
    {standardGeneric("reconstructGene")})

#' binary heatmap for standardized feature matrix
#' @export
#' @docType methods
#' @rdname binaryA-methods
#'
#' @description creates a binarized heatmap of the A matrix
#' in which the value is 1 if the value in Amean is greater than
#' threshold * Asd and 0 otherwise
#' @param object an object of type CogapsResult
#' @param threshold the number of standard deviations above zero
#' that an element of Amean must be to get a value of 1
#' @return plots a heatmap of the A Matrix
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' binMatrix <- binaryA(result, threshold=3)
setGeneric("binaryA", function(object, threshold=3)
    {standardGeneric("binaryA")})

#' plot of residuals
#' @export
#' @docType methods
#' @rdname plotResiduals-methods
#'
#' @description calculate residuals and produce heatmap
#' @param object an object of type CogapsResult
#' @param data original data matrix run through GAPS
#' @param uncertainty original standard deviation matrix run through GAPS
#' @return creates a residual plot
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' plotResiduals(result, SimpSim.D)
setGeneric("plotResiduals", function(object, data, uncertainty=NULL)
    {standardGeneric("plotResiduals")})

#' calculate gene set statistics
#' @export
#' @docType methods
#' @rdname calcCoGAPSStat-methods
#'
#' @description calculates the gene set statistics for each
#' column of A using a Z-score from the elements of the A matrix,
#' the input gene set, and permutation tests
#' @param object an object of type CogapsResult
#' @param GStoGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @return gene set statistics for each column of A
#' @examples
#' data('SimpSim')
#' result <- CoGAPS(SimpSim.D)
#' calcCoGAPSStat(result, GStoGenes=GSets, numPerm=500)
setGeneric("calcCoGAPSStat", function(object, GStoGenes, numPerm=500)
    {standardGeneric("calcCoGAPSStat")})

#' probability gene belongs in gene set
#' @export
#' @docType methods
#' @rdname calcGeneGSStat-methods
#'
#' @description calculates the probability that a gene
#' listed in a gene set behaves like other genes in the set within
#' the given data set
#' @param object an object of type CogapsResult
#' @param GSGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @param Pw weight on genes
#' @param nullGenes logical indicating gene adjustment
#' @return gene similiarity statistic
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' calcGeneGSStat(result, GSGenes=GSets[[1]], numPerm=500)
setGeneric("calcGeneGSStat", function(object, GStoGenes, numPerm,
Pw=rep(1,ncol(Amean)), nullGenes=FALSE)
    {standardGeneric("calcGeneGSStat")})

#' compute gene probability
#' @export
#' @docType methods
#' @rdname computeGeneGSProb-methods
#'
#' @description Computes the p-value for gene set membership using the CoGAPS-based
#' statistics developed in Fertig et al. (2012).  This statistic refines set
#' membership for each candidate gene in a set specified in \code{GSGenes} by
#' comparing the inferred activity of that gene to the average activity of the
#' set.
#' @param object an object of type CogapsResult
#' @param GSGenes data.frame or list with gene sets
#' @param Pw weight on genes
#' @param numPerm number of permutations for null
#' @param PwNull - logical indicating gene adjustment
#' @return A vector of length GSGenes containing the p-values of set membership
#' for each gene containined in the set specified in GSGenes.
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.D)
#' computeGeneGSProb(result, GSGenes=GSets[[1]], numPerm=500)
setGeneric("computeGeneGSProb", function(object, GStoGenes, numPerm=500,
Pw=rep(1,ncol(Amean)), PwNull=FALSE)
    {standardGeneric("computeGeneGSProb")})

################################## METHODS #####################################

setMethod("show", signature("CogapsResult"),
function(object)
{
    nFeatures <- nrow(object@featureLoadings)
    nSamples <- nrow(object@sampleFactors)
    nPatterns <- ncol(object@featureLoadings)

    print(paste("CogapsResult object with", nFeatures, "features and", nSamples,
        "samples"))
    print(paste(nPatterns, "patterns were learned"))
})

#' @export
#' @importFrom graphics plot
plot.CogapsResult <- function(object, ...)
{
    nSamples <- nrow(object@sampleFactors)
    nFactors <- ncol(object@sampleFactors)
    colors <- rainbow(nFactors)
    xlimits <- c(0, nSamples + 1)
    ylimits <- c(0, max(object@sampleFactors) * 1.1)

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude")

    for (i in 1:nFactors)
    {
        lines(x=1:nSamples, y=object@sampleFactors[,i], col=colors[i])
        points(x=1:nSamples, y=object@sampleFactors[,i], col=colors[i], pch=i)
    }

    legend("top", paste("Pattern", 1:nFactors, sep = ""), pch = 1:nFactors,
        lty=1, cex=0.8, col=colors, bty="y", ncol=5)
}

#' @rdname getMeanChiSq-methods
#' @aliases getMeanChiSq
setMethod("getMeanChiSq", signature(object="CogapsResult"),
function(object)
{
    object@metadata$meanChiSq
})

#' @rdname calcZ-methods
#' @aliases calcZ
setMethod("calcZ", signature(object="CogapsResult"),
function(object, which)
{
    if (which == "feature")
    {
        object@featureLoadings / object@featureStdDev
    }
    else if (which == "sample")
    {
        object@sampleFactors / object@sampleStdDev
    }
    else
    {
        stop("\"which\" must be either \"feature\" or \"sample\"")
    }
})

#' @rdname reconstructGene-methods
#' @aliases reconstructGene
setMethod("reconstructGene", signature(object="CogapsResult"),
function(object, genes)
{
    D <- object@featureLoadings %*% t(object@sampleFactors)
    if (!is.null(genes))
    {
        D <- D[genes,]
    }
    return(D)
})

#' @rdname binaryA-methods
#' @aliases binaryA
setMethod("binaryA", signature(object="CogapsResult"),
function(object, threshold)
{
    binA <- ifelse(calcZ(object) > threshold, 1, 0)

    heatmap.2(binA, Rowv = FALSE, Colv = FALSE, dendrogram="none",
        scale="none", col = brewer.pal(3,"Blues"), trace="none",
        density.info="none", cexCol=1.3, srtCol=45,
        lmat=rbind(c(0, 3), c(2,1), c(0,4) ),
        lwid=c(1,10), lhei=c(1, 4, 1.2 ),
        main="Heatmap of Standardized Feature Matrix")
    mtext(paste("(Threshold = ", threshold, ")"))
})

#' @rdname plotResiduals-methods
#' @aliases plotResiduals
setMethod("plotResiduals", signature(object="CogapsResult"),
function(object, data, uncertainty)
{
    data <- as.matrix(data)
    if (is.null(uncertainty))
        uncertainty <- pmax(0.1 * data, 0.1)
    uncertainty <- as.matrix(uncertainty)

    M <- reconstructGene(object)
    resid <- (data - M) / uncertainty
    
    scaledRdYlBu <- colorRampPalette(brewer.pal(9,"RdYlBu"))(100)
    heatmap.2(resid, Rowv = FALSE, Colv = FALSE, dendrogram="none",
        scale="none", col = scaledRdYlBu, trace="none", density.info="none",
        cexCol=1.33, srtCol=45, lmat=rbind(c(0, 3),c(2,1),c(0,4) ),
        lwid=c(1,10), lhei=c(1, 4, 1.2 ), main="Heatmap of Residuals")
})

#' patternMarkers
#'
#' @param threshold # the type of threshold to be used. The default "all" will distribute genes into pattern with the lowest ranking. The "cut" thresholding by the first gene to have a lower ranking, i.e. better fit to, a pattern.
#' @param lp a vector of weights for each pattern to be used for finding markers. If NA markers for each pattern of the A matrix will be used.
#' @param full logical indicating whether to return the ranks of each gene for each pattern
#' @return By default a non-overlapping list of genes associated with each \code{lp}. If \code{full=TRUE} a data.frame of
#' genes rankings with a column for each \code{lp} will also be returned.
#' @export

patternMarkers <- function(Amatrix=NA, scaledPmatrix=FALSE, Pmatrix=NA,
threshold="all", lp=NA, full=FALSE)

setMethod("patternMarkers", signature(object="CogapsResult"),
function(object, threshold="all", lp=NA, full=FALSE)

plotPatternMarkers

#' @rdname calcCoGAPSStat-methods
#' @aliases calcCoGAPSStat
setMethod("calcCoGAPSStat", signature(object="CogapsResult"),
function(object, GStoGenes, numPerm)
{
    # test for std dev of zero, possible mostly in simple simulations
    if (sum(object@featureStdDev==0) > 0)
    {
        object@featureStdDev[object@featureStdDev==0] <- 1e-6
    }

    # calculate Z scores
    zMatrix <- calcZ(object)

    # check input arguments
    if (!is(GStoGenes, "data.frame") && !is(GStoGenes, "list") && !is(GStoGenes,"GSA.genesets"))
    {
        stop("GStoGenes must be a data.frame,GSA.genesets, or list with format specified in the users manual.")
    }

    if (is(GStoGenes, "GSA.genesets"))
    {
        names(GStoGenes$genesets) <- GStoGenes$geneset.names
        GStoGenes <- GStoGenes$genesets
    }

    if (is(GStoGenes, "list"))
    {
        GStoGenesList <- GStoGenes
    }
    else
    {
        GStoGenesList <- list()
        for (i in 1:dim(GStoGenes)[2])
        {
            GStoGenesList[[as.character(colnames(GStoGenes)[i])]] <- as.character(unique(GStoGenes[,i]))
        }
    }

    # get dimensions
    numGS   <- length(names(GStoGenesList))
    numPatt <- dim(zMatrix)[2]
    numG    <- dim(zMatrix)[1]+0.9999

    # initialize matrices
    statsUp   <- matrix(ncol = numGS, nrow = numPatt)
    statsDown <- matrix(ncol = numGS, nrow = numPatt)
    actEst    <- matrix(ncol = numGS, nrow = numPatt)
    results   <- list()
    zPerm     <- matrix(ncol=numPerm,nrow=numPatt)

    # do permutation test
    for (gs in 1:numGS)
    {
        genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
        index <- rownames(zMatrix) %in% genes
        zValues <- zMatrix[index,1]
        numGenes <- length(zValues)
        label <- as.character(numGenes)

        if (!any(names(results)==label))
        {
            for (p in 1:numPatt)
            {
                for (j in 1:numPerm)
                {
                    temp <- floor(runif(numGenes,1,numG))
                    temp2 <- zMatrix[temp,p]
                    zPerm[p,j] <- mean(temp2)
                }
            }
            results[[label]] <- zPerm
        }
    }

    # get p-value
    for (p in 1:numPatt)
    {
        for (gs in 1:numGS)
        {
            genes <- GStoGenesList[[names(GStoGenesList)[gs]]]
            index <- rownames(zMatrix) %in% genes
            zValues <- zMatrix[index,p]
            zScore <- mean(zValues)

            numGenes <- length(zValues)
            label <- as.character(numGenes)

            permzValues <- results[[label]][p,]
            ordering <- order(permzValues)
            permzValues <- permzValues[ordering]

            statistic <- sum(zScore > permzValues)
            statUpReg <- 1 - statistic/length(permzValues)
            statsUp[p,gs] <- max(statUpReg, 1/numPerm)

            statistic <- sum(zScore < permzValues)
            statDownReg <- 1 - statistic/length(permzValues)
            statsDown[p,gs] <- max(statDownReg,1/numPerm)

            activity <- 1 - 2*max(statUpReg, 1/numPerm)
            actEst[p,gs] <- activity
        }
    }

    # format output
    colnames(statsUp) <- names(GStoGenesList)
    colnames(statsDown) <- names(GStoGenesList)
    colnames(actEst) <- names(GStoGenesList)

    rownames(statsUp) <- colnames(zMatrix)
    rownames(statsDown) <- colnames(zMatrix)
    rownames(actEst) <- colnames(zMatrix)

    results[['GSUpreg']] <- statsUp
    results[['GSDownreg']] <- statsDown
    results[['GSActEst']] <- actEst

    return(results)
})

#' @rdname calcGeneGSStat-methods
#' @aliases calcGeneGSStat
setMethod("calcGeneGSStat", signature(object="CogapsResult"),
function(object, GStoGenes, numPerm, Pw, nullGenes)
{
    gsStat <- calcCoGAPSStat(object, data.frame(GSGenes), numPerm=numPerm)
    gsStat <- gsStat$GSUpreg
    gsStat <- -log(gsStat)

    if (!all(is.na(Pw)))
    {
        if (length(Pw) != length(gsStat))
        {
            stop('Invalid weighting')
        }
        gsStat <- gsStat*Pw
    }

    if (nullGenes)
    {
        ZD <- object@featureLoadings[setdiff(row.names(object@featureLoadings), GSGenes),] /
            object@featureStdDev[setdiff(row.names(object@featureLoadings), GSGenes),]
    }
    else
    {
        ZD <- object@featureLoadings[GSGenes,]/object@featureStdDev[GSGenes,]
    }
    outStats <- apply(sweep(ZD,2,gsStat,FUN="*"),1,sum) / (sum(gsStat))
    outStats <- outStats / apply(ZD,1,sum)
    outStats[which(apply(ZD,1,sum) < 1e-6)] <- 0

    if (sum(gsStat) < 1e-6)
    {
        return(0)
    }
    return(outStats)
})

#' @rdname computeGeneGSProb-methods
#' @aliases computeGeneGSProb
setMethod("computeGeneGSProb", signature(object="CogapsResult"),
function(object, GStoGenes, numPerm, Pw, PwNull)
{
    geneGSStat <- calcGeneGSStat(object, Pw=Pw, GSGenes=GSGenes,
        numPerm=numPerm)

    if (PwNull)
    {
        permGSStat <- calcGeneGSStat(object, GSGenes=GSGenes, numPerm=numPerm,
            Pw=Pw, nullGenes=TRUE)
    }
    else
    {
        permGSStat <- calcGeneGSStat(object, GSGenes=GSGenes, numPerm=numPerm,
            nullGenes=TRUE)
    }

    finalStats <- sapply(GSGenes, function(x)
        length(which(permGSStat > geneGSStat[x])) / length(permGSStat))

    return(finalStats)
})

