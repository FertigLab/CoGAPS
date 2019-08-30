#' convert list output from c++ code to a CogapsResult object
#' @keywords internal
#'
#' @param returnList list from cogaps_cpp
#' @param allParams list of all parameters
#' @return CogapsResult object
#' @importFrom utils packageVersion
createCogapsResult <- function(returnList, allParams)
{
    res <- new("CogapsResult",
        Amean       = returnList$Amean,
        Asd         = returnList$Asd,
        Pmean       = returnList$Pmean,
        Psd         = returnList$Psd,
        meanChiSq   = returnList$meanChiSq,
        geneNames   = returnList$geneNames,
        sampleNames = returnList$sampleNames,
        diagnostics = append(returnList$diagnostics,
            list("params"=allParams$gaps, "version"=utils::packageVersion("CoGAPS")))
    )
    validObject(res)
    return(res)
}

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
#' @importFrom graphics plot legend lines points
#' @importFrom grDevices rainbow
plot.CogapsResult <- function(x, ...)
{
    nSamples <- nrow(x@sampleFactors)
    nFactors <- ncol(x@sampleFactors)
    colors <- rainbow(nFactors)
    xlimits <- c(0, nSamples + 1)
    ylimits <- c(0, max(x@sampleFactors) * 1.1)

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude",
        xlab="Samples")

    for (i in 1:nFactors)
    {
        lines(x=1:nSamples, y=x@sampleFactors[,i], col=colors[i])
        points(x=1:nSamples, y=x@sampleFactors[,i], col=colors[i], pch=i)
    }

    legend("top", paste("Pattern", 1:nFactors, sep = ""), pch = 1:nFactors,
        lty=1, cex=0.8, col=colors, bty="y", ncol=5)
}

#' @rdname getFeatureLoadings-methods
#' @aliases getFeatureLoadings
setMethod("getFeatureLoadings", signature(object="CogapsResult"),
function(object)
{
    object@featureLoadings
})

#' @rdname getSampleFactors-methods
#' @aliases getSampleFactors
setMethod("getSampleFactors", signature(object="CogapsResult"),
function(object)
{
    object@sampleFactors
})

#' @rdname getMeanChiSq-methods
#' @aliases getMeanChiSq
setMethod("getMeanChiSq", signature(object="CogapsResult"),
function(object)
{
    object@metadata$meanChiSq
})

#' @rdname getVersion-methods
#' @aliases getVersion
setMethod("getVersion", signature(object="CogapsResult"),
function(object)
{
    object@metadata$version
})

#' @rdname getOriginalParameters-methods
#' @aliases getOriginalParameters
setMethod("getOriginalParameters", signature(object="CogapsResult"),
function(object)
{
    object@metadata$params
})

#' @rdname getUnmatchedPatterns-methods
#' @aliases getUnmatchedPatterns
setMethod("getUnmatchedPatterns", signature(object="CogapsResult"),
function(object)
{
    if (!is.null(object@metadata$unmatchedPatterns))
        return(object@metadata$unmatchedPatterns)

    message("this result was not generated with a call to GWCoGAPS or scCoGAPS")
    return(NULL)
})

#' @rdname getClusteredPatterns-methods
#' @aliases getClusteredPatterns
setMethod("getClusteredPatterns", signature(object="CogapsResult"),
function(object)
{
    if (!is.null(object@metadata$clusteredPatterns))
        return(object@metadata$clusteredPatterns)

    message("this result was not generated with a call to GWCoGAPS or scCoGAPS")
    return(NULL)
})

#' @rdname getCorrelationToMeanPattern-methods
#' @aliases getCorrelationToMeanPattern
setMethod("getCorrelationToMeanPattern", signature(object="CogapsResult"),
function(object)
{
    if (!is.null(object@metadata$CorrToMeanPattern))
        return(object@metadata$CorrToMeanPattern)

    message("this result was not generated with a call to GWCoGAPS or scCoGAPS")
    return(NULL)
})

#' @rdname getSubsets-methods
#' @aliases getSubsets
setMethod("getSubsets", signature(object="CogapsResult"),
function(object)
{
    if (!is.null(object@metadata$subsets))
        return(object@metadata$subsets)

    message("this result was not generated with a call to GWCoGAPS or scCoGAPS")
    return(NULL)
})

#' @rdname calcZ-methods
#' @aliases calcZ
setMethod("calcZ", signature(object="CogapsResult"),
function(object, whichMatrix)
{

    if (whichMatrix=="featureLoadings")
    {
        if (sum(object@featureStdDev==0) > 0)
            stop("zeros detected in the standard deviation matrix")
        return(object@featureLoadings / object@featureStdDev)
    }
    else if (whichMatrix=="sampleFactors")
    {
        if (sum(object@sampleStdDev==0) > 0)
            stop("zeros detected in the standard deviation matrix")
        return(object@sampleFactors / object@sampleStdDev)
    }
    else
    {
        stop("whichMatrix must be either 'featureLoadings' or 'sampleFactors'")
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
#' @importFrom gplots heatmap.2
#' @importFrom graphics mtext
#' @importFrom RColorBrewer brewer.pal
setMethod("binaryA", signature(object="CogapsResult"),
function(object, threshold)
{
    binA <- ifelse(calcZ(object) > threshold, 1, 0)

    gplots::heatmap.2(binA, Rowv = FALSE, Colv = FALSE, dendrogram="none",
        scale="none", col = brewer.pal(3,"Blues"), trace="none",
        density.info="none", cexCol=1.3, srtCol=45,
        lmat=rbind(c(0, 3), c(2,1), c(0,4) ),
        lwid=c(1,10), lhei=c(1, 4, 1.2 ),
        main="Heatmap of Standardized Feature Matrix")
    mtext(paste("(Threshold = ", threshold, ")"))
})

#' @rdname plotResiduals-methods
#' @aliases plotResiduals
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
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
    gplots::heatmap.2(resid, Rowv = FALSE, Colv = FALSE, dendrogram="none",
        scale="none", col = scaledRdYlBu, trace="none", density.info="none",
        cexCol=1.33, srtCol=45, lmat=rbind(c(0, 3),c(2,1),c(0,4) ),
        lwid=c(1,10), lhei=c(1, 4, 1.2 ), main="Heatmap of Residuals")
})

patternMarkersCalculation <- function(mat, threshold, lp)
{
    nMeasurements <- nrow(mat)
    nPatterns <- ncol(mat)
    rowMax <- t(apply(mat, 1, function(x) x / max(x)))
    
    if (!is.na(lp))
    {
        if (length(lp) != nPatterns)
        {
            warning("lp length must equal the number of columns of the Amatrix")
        }
        sstat <- apply(rowMax, 1, function(x) sqrt(t(x-lp) %*% (x-lp)))
        ssranks <- rank(sstat)
        ssgenes.th <- names(sort(sstat,decreasing=FALSE,na.last=TRUE))
    }
    else
    {
        # determine which measurements are most associated with each pattern
        sstat <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=dimnames(mat))
        ssranks <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=dimnames(mat))
        ssgenes <- matrix(NA, nrow=nMeasurements, ncol=nPatterns, dimnames=NULL)
        for (i in 1:nPatterns)
        {
            lp <- rep(0,nPatterns)
            lp[i] <- 1
            sstat[,i] <- unlist(apply(rowMax, 1, function(x) sqrt(t(x-lp) %*% (x-lp))))
            ssranks[,i] <- rank(sstat[,i])
            ssgenes[,i] <- names(sort(sstat[,i], decreasing=FALSE, na.last=TRUE))
        }
        if (threshold=="cut")
        {
            geneThresh <- sapply(1:nPatterns, function(x) min(which(ssranks[ssgenes[,x],x] > apply(ssranks[ssgenes[,x],],1,min))))
            ssgenes.th <- sapply(1:nPatterns, function(x) ssgenes[1:geneThresh[x],x])
        }
        else if (threshold=="all")
        {
            pIndx <- unlist(apply(sstat, 1, which.min))
            gBYp <- list()
            for (i in sort(unique(pIndx)))
            {
                gBYp[[i]] <- sapply(strsplit(names(pIndx[pIndx==i]),"[.]"),function(x) x[[1]][1])
            }
            ssgenes.th <- lapply(1:max(sort(unique(pIndx))), function(x)
                ssgenes[which(ssgenes[,x] %in% gBYp[[x]]),x])
        }
        else
        {
            stop("invalid thresold argument")
        }
    }
    return(list("PatternMarkers"=ssgenes.th, "PatternRanks"=ssranks,
        "PatternMarkerScores"=sstat))       
}

#' @rdname patternMarkers-methods
#' @aliases patternMarkers
setMethod("patternMarkers", signature(object="CogapsResult"),
function(object, threshold, lp)
{
    patternMarkersCalculation(object@featureLoadings, threshold, lp)
})

#' @rdname amplitudeMarkers-methods
#' @aliases amplitudeMarkers
setMethod("amplitudeMarkers", signature(object="CogapsResult"),
function(object, threshold, lp)
{
    patternMarkersCalculation(object@sampleFactors, threshold, lp)
})

#' heatmap of original data clustered by pattern markers statistic
#' @export
#' @docType methods
#' @rdname plotPatternMarkers-methods
#'
#' @param object an object of type CogapsResult
#' @param data the original data as a matrix
#' @param patternPalette a vector indicating what color should be used
#' for each pattern
#' @param sampleNames names of the samples to use for labeling 
#' @param samplePalette  a vector indicating what color should be used
#' for each sample
#' @param colDenogram logical indicating whether to display sample denogram
#' @param heatmapCol pallelet giving color scheme for heatmap
#' @param scale character indicating if the values should be centered and scaled
#' in either the row direction or the column direction, or none. The
#' default is "row". 
#' @param ... additional graphical parameters to be passed to \code{heatmap.2}
#' @return heatmap of the \code{data} values for the \code{patternMarkers}
#' @seealso  \code{\link{heatmap.2}}
#' @importFrom gplots bluered
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust as.dist cor
plotPatternMarkers <- function(object, data, patternPalette, sampleNames,
samplePalette=NULL, heatmapCol=bluered, colDenogram=TRUE, scale="row", ...)
{
    if (is.null(samplePalette))
        samplePalette <- rep("black", ncol(data))

    ### coloring of genes
    patternCols <- rep(NA, length(unlist(patternMarkers)))
    names(patternCols) <- unlist(patternMarkers)
    for (i in 1:length(patternMarkers))
    {
        patternCols[patternMarkers[[i]]] <- patternPalette[i]
    }

    gplots::heatmap.2(as.matrix(data[unlist(patternMarkers),],...),
        symkey=FALSE,
        symm=FALSE, 
        scale=scale,
        col=heatmapCol,
        distfun=function(x) as.dist((1-cor(t(x)))/2),
        hclustfun=function(x) stats::hclust(x,method="average"),
        Rowv=FALSE,
        Colv=colDenogram,
        trace='none',
        RowSideColors = as.character(patternCols[unlist(patternMarkers)]),
        labCol= sampleNames,
        cexCol=0.8,
        ColSideColors=as.character(samplePalette),
        rowsep=cumsum(sapply(patternMarkers,length))
    )
}

#' @rdname calcCoGAPSStat-methods
#' @aliases calcCoGAPSStat
#' @importFrom methods is
setMethod("calcCoGAPSStat", signature(object="CogapsResult"),
function(object, sets, whichMatrix, numPerm, ...)
{
    ## do this for backwards compatibility
    if (!is.null(list(...)$GStoGenes))
        sets <- GStoGenes
    ## check input arguments
    if (!is(sets, "list"))
        stop("sets must be a list of either measurements or samples")
    ## do a permutation test
    zMatrix <- calcZ(object, whichMatrix)
    pvalUpReg <- sapply(sets,
        function(thisSet)
        {
            numValues <- length(thisSet)
            lessThanCount <- rep(0, ncol(zMatrix))
            actualZScore <- colMeans(zMatrix[rownames(zMatrix) %in% thisSet,])
            for (n in 1:numPerm)
            {
                permutedIndices <- sample(1:nrow(zMatrix), size=numValues, replace=FALSE)
                permutedZScore <- colMeans(zMatrix[permutedIndices,])
                lessThanCount <- lessThanCount + (actualZScore < permutedZScore)
            }
            return(lessThanCount / numPerm)
        })
    ## calculate other quantities of interest, return a list
    pvalDownReg <- 1 - pvalUpReg # this is techinically >=, but == should rarely happen
    activityEstimate <- 1 - 2 * pvalUpReg
    return(list(
        'twoSidedPValue'=pmin(pvalDownReg, pvalUpReg),
        'GSUpreg'=pvalUpReg,
        'GSDownreg'=pvalDownReg,
        'GSActEst'=activityEstimate
    ))
})

#' @rdname calcGeneGSStat-methods
#' @aliases calcGeneGSStat
setMethod("calcGeneGSStat", signature(object="CogapsResult"),
function(object, GStoGenes, numPerm, Pw, nullGenes)
{
    gsStat <- calcCoGAPSStat(object, data.frame(GStoGenes), numPerm=numPerm)
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
        ZD <- object@featureLoadings[setdiff(row.names(object@featureLoadings), GStoGenes),] /
            object@featureStdDev[setdiff(row.names(object@featureLoadings), GStoGenes),]
    }
    else
    {
        ZD <- object@featureLoadings[GStoGenes,]/object@featureStdDev[GStoGenes,]
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
    geneGSStat <- calcGeneGSStat(object, Pw=Pw, GStoGenes=GStoGenes,
        numPerm=numPerm)

    if (PwNull)
    {
        permGSStat <- calcGeneGSStat(object, GStoGenes=GStoGenes, numPerm=numPerm,
            Pw=Pw, nullGenes=TRUE)
    }
    else
    {
        permGSStat <- calcGeneGSStat(object, GStoGenes=GStoGenes, numPerm=numPerm,
            nullGenes=TRUE)
    }

    finalStats <- sapply(GStoGenes, function(x)
        length(which(permGSStat > geneGSStat[x])) / length(permGSStat))

    return(finalStats)
})


