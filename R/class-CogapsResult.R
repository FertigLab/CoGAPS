#' CogapsResult
#' @export
#'
#' @slot sampleStdDev std dev of the sampled P matrices
#' @slot featureStdDev std dev of the sampled A matrices
#' @description Contains all output from Cogaps run
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
setClass("CogapsResult", contains="LinearEmbeddingMatrix", slots=c(
    sampleStdDev = "ANY",   # Psd
    featureStdDev = "ANY"   # Asd
))

#' Constructor for CogapsResult
#' @param .Object CogapsResult object
#' @param Amean mean of sampled A matrices
#' @param Pmean mean of sampled P matrices
#' @param Asd std dev of sampled A matrices
#' @param Psd std dev of sampled P matrices
#' @param meanChiSq mean value of ChiSq statistic
#' @param geneNames names of genes in data
#' @param sampleNames names of samples in data
#' @param diagnostics assorted diagnostic reports from the run
#' @param ... initial values for slots
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
function(.Object, Amean, Pmean, Asd, Psd, meanChiSq, geneNames,
sampleNames, diagnostics=NULL, ...)
{
    if (is.null(geneNames))
        warning("no gene names given")
    if (is.null(sampleNames))
        warning("no sample names given")

    .Object@featureLoadings <- Amean
    .Object@sampleFactors <- Pmean
    .Object@featureStdDev <- Asd
    .Object@sampleStdDev <- Psd

    patternNames <- paste("Pattern", 1:ncol(Amean), sep="_")

    rownames(.Object@featureLoadings) <- geneNames
    colnames(.Object@featureLoadings) <- patternNames

    rownames(.Object@featureStdDev) <- geneNames
    colnames(.Object@featureStdDev) <- patternNames

    rownames(.Object@sampleFactors) <- sampleNames
    colnames(.Object@sampleFactors) <- patternNames

    rownames(.Object@sampleStdDev) <- sampleNames
    colnames(.Object@sampleStdDev) <- patternNames

    .Object@metadata[["meanChiSq"]] <- meanChiSq
    .Object@metadata <- append(.Object@metadata, diagnostics)

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
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getMeanChiSq(result)
setGeneric("getMeanChiSq", function(object)
    {standardGeneric("getMeanChiSq")})

#' return version number used to generate this result
#' @export
#' @docType methods
#' @rdname getVersion-methods
#'
#' @param object an object of type CogapsResult
#' @return version number
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getVersion(result)
setGeneric("getVersion", function(object)
    {standardGeneric("getVersion")})

#' return original parameters used to generate this result
#' @export
#' @docType methods
#' @rdname getOriginalParameters-methods
#'
#' @param object an object of type CogapsResult
#' @return CogapsParams object
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getOriginalParameters(result)
setGeneric("getOriginalParameters", function(object)
    {standardGeneric("getOriginalParameters")})

#' return unmatched patterns from each subset
#' @export
#' @docType methods
#' @rdname getUnmatchedPatterns-methods
#'
#' @param object an object of type CogapsResult
#' @return CogapsParams object
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getUnmatchedPatterns(result)
setGeneric("getUnmatchedPatterns", function(object)
    {standardGeneric("getUnmatchedPatterns")})

#' return clustered patterns from set of all patterns across all subsets
#' @export
#' @docType methods
#' @rdname getClusteredPatterns-methods
#'
#' @param object an object of type CogapsResult
#' @return CogapsParams object
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getClusteredPatterns(result)
setGeneric("getClusteredPatterns", function(object)
    {standardGeneric("getClusteredPatterns")})

#' return correlation between each pattern and the cluster mean
#' @export
#' @docType methods
#' @rdname getCorrelationToMeanPattern-methods
#'
#' @param object an object of type CogapsResult
#' @return CogapsParams object
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getCorrelationToMeanPattern(result)
setGeneric("getCorrelationToMeanPattern", function(object)
    {standardGeneric("getCorrelationToMeanPattern")})

#' return the names of the genes (samples) in each subset
#' @export
#' @docType methods
#' @rdname getSubsets-methods
#'
#' @param object an object of type CogapsResult
#' @return CogapsParams object
#' @examples
#' data(SimpSim)
#' result <- CoGAPS(SimpSim.data)
#' getSubsets(result)
setGeneric("getSubsets", function(object)
    {standardGeneric("getSubsets")})

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
#' result <- CoGAPS(SimpSim.data)
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
#' result <- CoGAPS(SimpSim.data)
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
#' result <- CoGAPS(SimpSim.data)
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
#' result <- CoGAPS(SimpSim.data)
#' plotResiduals(result, SimpSim.data)
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
#' @param GStoGenes data.frame or list with gene sets
#' @param numPerm number of permutations for null
#' @param Pw weight on genes
#' @param nullGenes logical indicating gene adjustment
#' @return gene similiarity statistic
setGeneric("calcGeneGSStat", function(object, GStoGenes, numPerm,
Pw=rep(1,ncol(object@featureLoadings)), nullGenes=FALSE)
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
#' @param GStoGenes data.frame or list with gene sets
#' @param Pw weight on genes
#' @param numPerm number of permutations for null
#' @param PwNull - logical indicating gene adjustment
#' @return A vector of length GSGenes containing the p-values of set membership
#' for each gene containined in the set specified in GSGenes.
setGeneric("computeGeneGSProb", function(object, GStoGenes, numPerm=500,
Pw=rep(1,ncol(object@featureLoadings)), PwNull=FALSE)
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
setGeneric("patternMarkers", function(object, threshold="all", lp=NA)
    {standardGeneric("patternMarkers")})

