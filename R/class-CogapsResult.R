#' CogapsResult
#' @export
#' @example 
#' # Return output from CogapsResult
#' CogapsResult <- CoGAPS(GIST.D, CogapsParams)
#' @description Contains all output from Cogaps run
setClass("CogapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix",
    seed = "numeric",
    meanChiSq = "numeric",
    diagnostics = "list"
))

#' Constructor for CogapsResult
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
function(.Object, ...)
{
    .Object <- callNextMethod(.Object, ...)
    .Object
})

setMethod("show", signature("CogapsResult"),
function(object)
{
    nGenes <- nrow(object@Amean)
    nPatterns <- ncol(object@Amean)
    nSamples <- ncol(object@Pmean)

    print(paste("CogapsResult object with", nGenes, "genes and", nSamples, "samples"))
    print(paste(nPatterns, "patterns were learned"))
})

setMethod('plot', signature('CoGAPSResult'),
function(object)
{
    colors <- rainbow(nrow(object@Pmean))
    xlimits <- c(0, ncol(object@Pmean) + 1)
    ylimits <- c(0, (max(object@Pmean) * 1.05))

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude")

    for (i in 1:nrow(object@Pmean))
    {
        lines(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i])
        points(x=1:ncol(object@Pmean), y=object@Pmean[i,], col=colors[i], pch=i)
    }

    legend("bottom", paste("Pattern", 1:nrow(object@Pmean), sep = ""),
    pch = 1:nrow(object@Pmean), lty=1, cex=0.8, col=colors, bty="y", ncol=5)
})

setMethod("MergeResultsWithSCE", signature("CogapsResult", "SingleCellExperiment"),
function(result, SCE)
{
    SCE@reducedDims <- SimpleList(Amean=result@Amean, Pmean=result@Pmean)
    return(SCE)
})

