#' CogapsResult
#' @export
#'
#' @description Contains all output from Cogaps run
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
setClass("CogapsResult", contains="LinearEmbeddingMatrix", slots=c(
    sampleStdDev = "ANY",   # Psd transpose
    featureStdDev = "ANY"   # Asd
))

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

#' compute pattern markers statistic
#' @export
#' @docType methods
#' @rdname patternMarkers-methods
#'
#' @description calculate the most associated pattern for each gene
#' @param object an object of type CogapsResult
#' @param threshold the type of threshold to be used. The default "all" will
#' distribute genes into pattern with the lowest ranking. The "cut" thresholds
#' by the first gene to have a lower ranking, i.e. better fit to, a pattern.
#' @param lp a vector of weights for each pattern to be used for finding
#' markers. If NA markers for each pattern of the A matrix will be used.
#' @return By default a non-overlapping list of genes associated with each
#' \code{lp}. If \code{full=TRUE} a data.frame of genes rankings with a column
#' for each \code{lp} will also be returned.
setGeneric("patternMarkers", function(object, threshold="all", lp=NA, full=FALSE)
    {standardGeneric("patternMarkers")})

