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

    # label snapshots
    if (allParams$nSnapshots > 0)
    {
        freq <- floor(allParams$gaps@nIterations / allParams$nSnapshots)
        labels <- seq(freq, allParams$gaps@nIterations, freq)
        if (allParams$snapshotPhase == 'equilibration' | allParams$snapshotPhase == 'all')
        {
            for (n in 1:allParams$nSnapshots)
            {
                dimnames(res@metadata$equilibrationSnapshotsA[[n]]) <- dimnames(res@featureLoadings)
                dimnames(res@metadata$equilibrationSnapshotsP[[n]]) <- dimnames(res@sampleFactors)
            }
            names(res@metadata$equilibrationSnapshotsA) <- labels
            names(res@metadata$equilibrationSnapshotsP) <- labels
        }
        if (allParams$snapshotPhase == 'sampling' | allParams$snapshotPhase == 'all')
        {
            for (n in 1:allParams$nSnapshots)
            {
                dimnames(res@metadata$samplingSnapshotsA[[n]]) <- dimnames(res@featureLoadings)
                dimnames(res@metadata$samplingSnapshotsP[[n]]) <- dimnames(res@sampleFactors)
            }
            names(res@metadata$samplingSnapshotsA) <- labels
            names(res@metadata$samplingSnapshotsP) <- labels
        }
    }        
    
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
#' @importFrom graphics plot legend lines points axis par text
#' @importFrom grDevices rainbow
plot.CogapsResult <- function(x, groups=NULL, ...)
{
  
  if(!is.null(groups)) {
    grpset <- unique(groups)
    groups <- table(groups)
    ngroups <- length(grpset)
    nFactors <- ncol(x@sampleFactors)
    colors <- rainbow(nFactors)
    xlimits <- c(0, ngroups + 1)
    ylimits <- c(0, max(x@sampleFactors) * 1.1)
    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude", xlab="", axes=FALSE) 
    axis(1, labels=FALSE, at=1:length(grpset))
    text(x = seq(1, length(grpset), by=1), par("usr")[3]-0.15, labels = grpset, srt = 60, pos = 1, xpd = TRUE)
    axis(2)
    
    for (i in 1:nFactors)
    {
      lines(x=1:length(groups), y=x@sampleFactors[groups,i], col=colors[i])
      points(x=1:length(groups), y=x@sampleFactors[groups,i], col=colors[i], pch=i)
    }
    
    legend("top", paste("Pattern", 1:nFactors, sep = ""), pch = 1:nFactors,
           lty=1, cex=0.8, col=colors, bty="y", ncol=5)
  }
  else{
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
}

#' @rdname getFeatureLoadings-methods
#' @aliases getFeatureLoadings
setMethod("getFeatureLoadings", signature(object="CogapsResult"),
function(object)
{
    object@featureLoadings
})



#' @rdname getAmplitudeMatrix-methods
#' @aliases getAmplitudeMatrix
setMethod("getAmplitudeMatrix", signature(object="CogapsResult"),
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

#' @rdname getPatternMatrix-methods
#' @aliases getPatternMatrix
setMethod("getPatternMatrix", signature(object="CogapsResult"),
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
    if (!(whichMatrix %in% c("featureLoadings", "sampleFactors")))
        stop("whichMatrix must be either 'featureLoadings' or 'sampleFactors'")
    mean <- if (whichMatrix == "sampleFactors") object@sampleFactors else object@featureLoadings
    stddev <- if (whichMatrix == "sampleFactors") object@factorStdDev else object@loadingStdDev
    if (sum(stddev == 0) > 0)
        warning("zeros detected in the standard deviation matrix")
    stddev[stddev == 0] <- 1e-6
    return(mean / stddev)    
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

unitVector <- function(n, length)
{
    vec <- rep(0, length)
    vec[n] <- 1
    return(vec)
}

#' @rdname getPatternGeneSet-methods
#' @import dplyr
#' @import fgsea
#' @importFrom forcats fct_reorder
#' @aliases getPatternGeneSet
setMethod("getPatternGeneSet", signature(object="CogapsResult", gene.sets="list", method = "character"),
function(object, gene.sets, method = c("enrichment", "overrepresentation"), ...)
{
  method <- match.arg(method)
  A <- object@featureLoadings
  patternNames <- colnames(A)
  features <- rownames(A)
  # enrichment method
  if(method == "enrichment") {
    d <- lapply(
      patternNames, FUN = function(p) {
        amp <- A[,p]
        names(amp) <- features
        result <- fgsea::fgsea(pathways = gene.sets, stats = amp, scoreType = "pos")
        result$leadingEdge <- vapply(result$leadingEdge, FUN = toString, FUN.VALUE = character(1))
        result$neg.log.padj <- (-10) * log10(result$padj)
        result$gene.set <- names(gene.sets)
        result <- dplyr::mutate(result,gene.set=forcats::fct_reorder(gene.set, - padj))
        return(result)
      }
    )
  }
  
  # overrepresentation method
  if(method == "overrepresentation") {
    patternMarkerResults <- patternMarkers(object, ...)
    names(patternMarkerResults$PatternMarkers) <- patternNames
    
    PMlist <- patternMarkerResults$PatternMarkers
    d <- lapply(
      patternNames, FUN = function(p) {
        result <- fgsea::fora(
          pathways = gene.sets, genes = PMlist[[p]],
          universe = features,
          maxSize=2038)
        result[["k/K"]] <- result$overlap / result$size
        result$neg.log.padj <- (-10) * log10(result$padj)
        result$gene.set <- names(gene.sets)
        result <- dplyr::mutate(result,gene.set=forcats::fct_reorder(gene.set, - padj))
        return(result)
      }
    )
  }
  return(d)
})

#' @rdname plotPatternGeneSet-methods
#' @importFrom dplyr relocate
#' @importFrom ggplot2 ggplot geom_col ylab coord_flip theme_minimal ggtitle geom_hline geom_text theme aes ylim
#' @importFrom graphics plot legend lines points
#' @aliases plotPatternGeneSet
setMethod("plotPatternGeneSet", signature(patterngeneset = "list", whichpattern="numeric", padj_threshold = "numeric"),
function(patterngeneset, whichpattern=1, padj_threshold = 0.05)
{
  # should be able to pass in info about just one pattern, or to pass all info and specify which pattern
  if (length(patterngeneset) > 1){
    gs_df <- patterngeneset[[whichpattern]]
  } else {
    gs_df <- patterngeneset[1]
  }
  # check the test type conducted to inform the title of the resulting plot
  if("k/K" %in% colnames(gs_df)) {
    method_name <- "Overrepresented"
  } else {
    method_name <- "Enriched"
  }
  
  # for plotting, limit the results to significant over-representation
  gs_df <- gs_df[gs_df$padj < padj_threshold,]
  neg.log.hline <- -10*log10(0.05)
  
  axis_max <- ceiling(max(gs_df$neg.log.padj)) + (max(gs_df$neg.log.padj)/4)
  
  #plot and save the waterfall plot of ORA p-values
  plot <- ggplot(gs_df, aes(y = neg.log.padj, x = gene.set, fill = gene.set)) +
    geom_col() +
    ylab("-10*log10(FDR-adjusted p-value)") + 
    coord_flip() +
    theme_minimal() +
    ggtitle(paste0(method_name, " gene sets in Pattern_", whichpattern)) +
    geom_text(aes(label=format(signif(padj, 4))), hjust = -.04) +
    theme(legend.position = "none") +
    ylim(0, axis_max)
  
  if(neg.log.hline <= axis_max) {
    plot <- plot +
      geom_hline(yintercept=neg.log.hline, linetype="dotted")
  }
  
  return(plot)
})


#' @rdname patternMarkers-methods
#' @aliases patternMarkers
setMethod("patternMarkers", signature(object="CogapsResult"),
function(object, threshold, lp, axis)
{
    ## check inputs to the function
    if (!(threshold %in% c("cut", "all")))
        stop("threshold must be either 'cut' or 'all'")
    if (!is.na(lp) & length(lp) != ncol(object@featureLoadings))
        stop("lp length must equal the number of patterns")
    if (!(axis %in% 1:2))
        stop("axis must be either 1 or 2")
    ## need to scale each row of the matrix of interest so that the maximum is 1
    resultMatrix <- if (axis == 1) object@featureLoadings else object@sampleFactors
    normedMatrix <- t(apply(resultMatrix, 1, function(row) row / max(row)))
    ## handle the case where a custom linear combination of patterns was passed in
    if (!is.na(lp))
    {
        markerScores <- apply(normedMatrix, 1, function(row) sqrt(sum((row-lp)^2)))
        markersByPattern <- names(sort(markerScores, decreasing=FALSE, na.last=TRUE))
        return(list(
            "PatternMarkers"=markersByPattern,
            "PatternMarkerRanks"=rank(markerScores),
            "PatternMarkerScores"=markerScores
        ))
    }
    ## default pattern marker calculation, each pattern has unit weight
    markerScores <- sapply(1:ncol(normedMatrix), function(patternIndex)
        apply(normedMatrix, 1, function(row)
        {
            lp <- unitVector(patternIndex, ncol(normedMatrix))
            return(sqrt(sum((row-lp)^2)))
        })
    )
    markerRanks <- apply(markerScores, 2, rank)
    colnames(markerScores) <- colnames(markerRanks) <- colnames(normedMatrix)
    ## keep only a subset of markers for each pattern depending on the type of threshold
    if (threshold == "cut") # all markers which achieve minimum rank
    {
      simplicityGENES <- function(As, Ps) {
        # rescale p's to have max 1
        pscale <- apply(Ps,1,max)
        
        # rescale A in accordance with p's having max 1
        As <- sweep(As, 2, pscale, FUN="*")
        
        # find the A with the highest magnitude
        Arowmax <- t(apply(As, 1, function(x) x/max(x)))
        pmax <- apply(As, 1, max)
        
        # determine which genes are most associated with each pattern
        ssl <- matrix(NA, nrow=nrow(As), ncol=ncol(As),
                      dimnames=dimnames(As))
        for (i in 1:ncol(As)) {
          lp <- rep(0, ncol(As))
          lp[i] <- 1
          ssl.stat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
          ssl[order(ssl.stat),i] <- 1:length(ssl.stat)
        }
        
        return(ssl)
        
      }
      simGenes <- simplicityGENES(As=object@featureLoadings,
                                  Ps=object@sampleFactors)
      
      patternMarkers <- list()
      
      nP <- ncol(simGenes)
      
      for (i in 1:nP) {
        
        sortSim <- names(sort(simGenes[,i],decreasing=F))
        
        geneThresh <- min(which(simGenes[sortSim,i] > 
                                  apply(simGenes[sortSim,],1,min)))
        
        markerGenes <- sortSim[1:geneThresh]
        
        markerGenes <- unique(markerGenes)
        
        patternMarkers[[i]] <- markerGenes
        
        markersByPattern <- patternMarkers
        
      }
    }
    else if (threshold == "all") # only the markers with the lowest scores
    {
        thresholdTest <- apply(markerScores, 1, which.min)
        patternsByMarker <- rep(0, length(markerScores))
        i <- 0
        for (gene in thresholdTest) {
          if (length(names(gene))){
            patternsByMarker[i] <- names(gene)
          }          
          i <- i + 1
        }

        markersByPattern <- sapply(colnames(markerScores), USE.NAMES=TRUE, simplify=FALSE,
            function(pattern) rownames(markerScores)[which(patternsByMarker==pattern)])
        
    }
    return(list(
        "PatternMarkers"=markersByPattern,
        "PatternMarkerRanks"=markerRanks,
        "PatternMarkerScores"=markerScores
    ))
})

#' @rdname calcCoGAPSStat-methods
#' @aliases calcCoGAPSStat
#' @importFrom methods is
setMethod("calcCoGAPSStat", signature(object="CogapsResult"),
function(object, sets, whichMatrix, numPerm, ...)
{
    ## do this for backwards compatibility
    if (!is.null(list(...)$GStoGenes))
        sets <- list(...)$GStoGenes
    ## check input arguments
    if (!is(sets, "list"))
        stop("sets must be a list of either measurements or samples")
    ## do a permutation test
    zMatrix <- calcZ(object, whichMatrix)
    pvalUpReg <- sapply(sets, function(thisSet)
    {
        lessThanCount <- rep(0, ncol(zMatrix))
        actualZScore <- colMeans(zMatrix[rownames(zMatrix) %in% thisSet,])
        for (n in 1:numPerm)
        {
            permutedIndices <- sample(1:nrow(zMatrix), size=length(thisSet), replace=FALSE)
            permutedZScore <- colMeans(zMatrix[permutedIndices,])
            lessThanCount <- lessThanCount + (actualZScore < permutedZScore)
        }
        return(lessThanCount / numPerm)
    })
    ## calculate other quantities of interest, return a list
    pvalDownReg <- 1 - pvalUpReg # this is techinically >=, but == should rarely happen
    activityEstimate <- 1 - 2 * pvalUpReg
    return(list(
        'twoSidedPValue'=pmax(pmin(pvalDownReg, pvalUpReg), 1 / numPerm),
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
            object@loadingStdDev[setdiff(row.names(object@featureLoadings), GStoGenes),]
    }
    else
    {
        ZD <- object@featureLoadings[GStoGenes,]/object@loadingStdDev[GStoGenes,]
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

#' @rdname MANOVA-methods
#' @aliases MANOVA
#' @importFrom stats manova
setMethod("MANOVA", signature(interestedVariables = "matrix", object = "CogapsResult"), 
function(interestedVariables, object){
  interestedVariables <- cbind(unclass(factor(interestedVariables[,1])), unclass(factor(interestedVariables[,2])))
  
  pat <- as.data.frame(object@sampleFactors)
  npat <- ncol(pat)
  pattern_names = colnames(pat)
  
  fits <- list()
  
  for (pattern in pattern_names) {
    print(pattern)
    pattern_column <- unlist(pat[,pattern])
    
    fit <- manova(interestedVariables ~ pattern_column)
    fits[[length(fits)+1]] <- fit
    print(summary(fit))
  }
  names(fits)=pattern_names
  return(fits)
})

#' @rdname toCSV-methods
#' @aliases toCSV
#' @importFrom utils write.csv
setMethod("toCSV", signature(object="CogapsResult", save_location="character"),
function(object, save_location)
{
  write.csv(object@featureLoadings, file = paste0(save_location, "/featureLoadings.csv"), row.names=FALSE)
  write.csv(object@sampleFactors, file= paste0(save_location, "/sampleFactors.csv"), row.names = FALSE)
  
  write.csv(object@loadingStdDev, file = paste0(save_location, "/loadingStdDev.csv"), row.names = FALSE)
  write.csv(object@factorStdDev, file = paste0(save_location, "/factorStdDev.csv"), row.names = FALSE)
  
  write.csv(rownames(object@featureLoadings), file = paste0(save_location, "/geneNames.csv"), row.names = FALSE)
  write.csv(rownames(object@sampleFactors), file = paste0(save_location, "/sampleNames.csv"), row.names = FALSE)
  
  write.csv(object@metadata$meanChiSq, file=paste0(save_location, "/meanChiSq.csv"), row.names = FALSE)
  write.csv(object@metadata$chisq, file = paste0(save_location, "/chisqHistory.csv"), row.names = FALSE)
  
  write.csv(object@metadata$version, file = paste0(save_location, "/version.csv"), row.names = FALSE)
  write.csv(object@metadata$atomsA, file = paste0(save_location, "/atomsA.csv"), row.names = FALSE)
  write.csv(object@metadata$atomsP, file = paste0(save_location, "/atomsP.csv"), row.names = FALSE)
})

#' @rdname fromCSV-methods
#' @aliases fromCSV
#' @importFrom utils read.csv
setMethod("fromCSV", signature(save_location="character"),
function(save_location)
{
  featureLoadings <- read.csv(file = paste0(save_location, "/featureLoadings.csv"))
  sampleFactors <- read.csv(file = paste0(save_location, "/sampleFactors.csv"))
  
  loadingStdDev <- read.csv(file = paste0(save_location, "/loadingStdDev.csv"))
  factorStdDev <- read.csv(file = paste0(save_location, "/factorStdDev.csv"))
  
  geneNames <- read.csv(file=paste0(save_location, "/geneNames.csv"))$x
  sampleNames <- read.csv(file = paste0(save_location, "/sampleNames.csv"))$x
  
  meanChiSq <- read.csv(file=paste0(save_location, "/meanChiSq.csv"))$x
  chisqHistory <- read.csv(file=paste0(save_location, "/chisqHistory.csv"))$x
  
  version <- read.csv(file=paste0(save_location, "/version.csv"))$x
  atomsA <- read.csv(file=paste0(save_location, "/atomsA.csv"))$x
  atomsP <- read.csv(file=paste0(save_location, "/atomsP.csv"))$x
  
  res <- new("CogapsResult",
      Amean       = featureLoadings,
      Asd         = loadingStdDev,
      Pmean       = sampleFactors,
      Psd         = factorStdDev,
      meanChiSq   = meanChiSq,
      geneNames   = geneNames,
      sampleNames = sampleNames
  )
  
  res@metadata$chisq = chisqHistory
  res@metadata$version = version
  res@metadata$atomsA = atomsA
  res@metadata$atomsP = atomsP
  
  return(res)
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
#' @param patternMarkers pattern markers to be plotted, as generated by the 
#' patternMarkers function
#' @param sampleNames names of the samples to use for labeling 
#' @param samplePalette  a vector indicating what color should be used
#' for each sample
#' @param colDendrogram logical indicating whether to display sample dendrogram
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
plotPatternMarkers <- function(object, data, patternMarkers, patternPalette, sampleNames,
samplePalette=NULL, heatmapCol=bluered, colDendrogram=TRUE, scale="row", ...)
{
    if (is.null(samplePalette))
        samplePalette <- rep("black", ncol(data))

    ### coloring of genes
    patternCols <- rep(NA, length(unlist(patternMarkers)))
    names(patternCols) <- unlist(patternMarkers)
    for (i in 1:length(patternMarkers))
    {
        patternCols[unlist(patternMarkers$PatternMarkers[[i]])] <- patternPalette[i]
    }

    gplots::heatmap.2(as.matrix(data[unlist(patternMarkers$PatternMarkers),],),
        symkey=FALSE,
        symm=FALSE, 
        scale=scale,
        col=heatmapCol,
        distfun=function(x) as.dist((1-cor(t(x)))/2),
        hclustfun=function(x) stats::hclust(x,method="average"),
        Rowv=FALSE,
        Colv=colDendrogram,
        trace='none',
        RowSideColors = as.character(patternCols[unlist(patternMarkers$PatternMarkers)]),
        labCol= sampleNames,
        cexCol=0.8,
        ColSideColors=as.character(samplePalette),
        rowsep=cumsum(sapply(patternMarkers,length))
    )
}
