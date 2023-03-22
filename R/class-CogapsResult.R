#' CogapsResult
#' @export
#'
#' @slot factorStdDev std dev of the sampled P matrices
#' @slot loadingStdDev std dev of the sampled A matrices
#' @description Contains all output from Cogaps run
#' @importClassesFrom S4Vectors Annotated
#' @importClassesFrom SingleCellExperiment LinearEmbeddingMatrix
setClass("CogapsResult", contains="LinearEmbeddingMatrix", slots=c(
    factorStdDev = "ANY",   # Psd
    loadingStdDev = "ANY"   # Asd
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
#' @importFrom S4Vectors make_zero_col_DFrame
setMethod("initialize", "CogapsResult",
function(.Object, Amean, Pmean, Asd, Psd, meanChiSq, geneNames,
sampleNames, diagnostics=NULL, ...)
{
    if (is.null(geneNames))
        stop("no gene names given")
    if (is.null(sampleNames))
        stop("no sample names given")

    .Object@featureLoadings <- Amean
    .Object@sampleFactors <- Pmean
    .Object@loadingStdDev <- Asd
    .Object@factorStdDev <- Psd

    patternNames <- paste("Pattern", 1:ncol(Amean), sep="_")

    if (length(geneNames) != nrow(.Object@featureLoadings))
        stop("number of gene names doesn't match data size, ",
            length(geneNames), " != ", nrow(.Object@featureLoadings))
    if (length(sampleNames) != nrow(.Object@sampleFactors))
        stop("number of sample names doesn't match data size")

    rownames(.Object@featureLoadings) <- geneNames
    colnames(.Object@featureLoadings) <- patternNames

    rownames(.Object@loadingStdDev) <- geneNames
    colnames(.Object@loadingStdDev) <- patternNames

    rownames(.Object@sampleFactors) <- sampleNames
    colnames(.Object@sampleFactors) <- patternNames

    rownames(.Object@factorStdDev) <- sampleNames
    colnames(.Object@factorStdDev) <- patternNames

    .Object@metadata[["meanChiSq"]] <- meanChiSq
    .Object@metadata <- append(.Object@metadata, diagnostics)

    .Object@factorData <- make_zero_col_DFrame(ncol(.Object@sampleFactors))

    .Object <- callNextMethod(.Object, ...)
    .Object
})

setValidity("CogapsResult",
    function(object)
    {
        if (any(is.na(object@featureLoadings)) | any(object@featureLoadings == Inf) | any(object@featureLoadings == -Inf))
            "NA/Inf values in feature matrix"
        if (any(is.na(object@sampleFactors)) | any(object@sampleFactors == Inf) | any(object@sampleFactors == -Inf))
            "NA/Inf values in sample matrix"
        # TODO: fix this check for semisupervised case
    #  if (sum(object@featureLoadings < 0) > 0 | sum(object@loadingStdDev < 0) > 0)
     #       "negative values in feature Matrix"
      #  if (sum(object@sampleFactors < 0) > 0 | sum(object@factorStdDev < 0) > 0)
       #     "negative values in sample Matrix"
    }
)    

################################### GENERICS ###################################

#' return featureLoadings matrix from CogapsResult object
#' @export
#' @docType methods
#' @rdname getFeatureLoadings-methods
#'
#' @param object an object of type CogapsResult
#' @return featureLoadings matrix
#' @examples
#' data(GIST)
#' fLoadings <- getFeatureLoadings(GIST.result)
setGeneric("getFeatureLoadings", function(object) standardGeneric("getFeatureLoadings"))

#' return Amplitude  matrix from CogapsResult object
#' @export
#' @docType methods
#' @rdname getAmplitudeMatrix-methods
#'
#' @param object an object of type CogapsResult
#' @return amplitude matrix
#' @examples
#' data(GIST)
#' amplitudeMatrix <- getAmplitudeMatrix(GIST.result)
setGeneric("getAmplitudeMatrix", function(object)
    {standardGeneric("getAmplitudeMatrix")})

#' return sampleFactors matrix from CogapsResult object
#' @export
#' @docType methods
#' @rdname getSampleFactors-methods
#'
#' @param object an object of type CogapsResult
#' @return sampleFactors matrix
#' @examples
#' data(GIST)
#' sFactors <- getSampleFactors(GIST.result)
setGeneric("getSampleFactors", function(object)
    {standardGeneric("getSampleFactors")})



#' return pattern matrix from CogapsResult object
#' @export
#' @docType methods
#' @rdname getPatternMatrix-methods
#'
#' @param object an object of type CogapsResult
#' @return pattern matrix
#' @examples
#' data(GIST)
#' patternMatrix <- getPatternMatrix(GIST.result)
setGeneric("getPatternMatrix", function(object)
    {standardGeneric("getPatternMatrix")})


#' return chi-sq of final matrices
#' @export
#' @docType methods
#' @rdname getMeanChiSq-methods
#'
#' @param object an object of type CogapsResult
#' @return chi-sq error
#' @examples
#' data(GIST)
#' getMeanChiSq(GIST.result)
setGeneric("getMeanChiSq", function(object)
    {standardGeneric("getMeanChiSq")})

#' generate statistics associating patterns with MSigDB hallmark gene sets
#' @export
#' @docType methods
#' @rdname getPatternHallmarks-methods
#' @aliases getPatternHallmarks
#' @param object an object of type CogapsResult
#' @return dataframe of hallmark info
setGeneric("getPatternHallmarks", function(object) standardGeneric("getPatternHallmarks"))

#' generate a barchart of most significant hallmark sets for a pattern
#' @export
#' @docType methods
#' @rdname plotPatternHallmarks-methods
#' @aliases plotPatternHallmarks
#' @param object an object of type CogapsResult
#' @param patternhallmarks output from getPatternHallmarks
#' @param whichpattern which pattern to generate bar chart for
#' @return image object of barchart
setGeneric("plotPatternHallmarks", function(object, patternhallmarks, whichpattern=1) standardGeneric("plotPatternHallmarks"))

#' return version number used to generate this result
#' @export
#' @docType methods
#' @rdname getVersion-methods
#'
#' @param object an object of type CogapsResult
#' @return version number
#' @examples
#' data(GIST)
#' getVersion(GIST.result)
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
#' data(GIST)
#' params <- getOriginalParameters(GIST.result)
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
#' data(GIST)
#' unmatchedPatterns <- getUnmatchedPatterns(GIST.result)
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
#' data(GIST)
#' clusteredPatterns <- getClusteredPatterns(GIST.result)
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
#' data(GIST)
#' corrToMeanPattern <- getCorrelationToMeanPattern(GIST.result)
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
#' data(GIST)
#' subsets <- getSubsets(GIST.result)
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
#' @param whichMatrix either "featureLoadings" or "sampleFactors" indicating which
#' matrix to calculate the z-score for
#' @return matrix of z-scores
#' @examples
#' data(GIST)
#' featureZScore <- calcZ(GIST.result, "featureLoadings")
setGeneric("calcZ", function(object, whichMatrix)
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
#' data(GIST)
#' estimatedD <- reconstructGene(GIST.result)
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
#' data(GIST)
#' # to expensive to call since it plots
#' # binaryA(GIST.result, threshold=3)
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
#' data(GIST)
#' # to expensive to call since it plots
#' # plotResiduals(GIST.result, GIST.matrix)
setGeneric("plotResiduals", function(object, data, uncertainty=NULL)
    {standardGeneric("plotResiduals")})

#' calculate statistic on sets of measurements (genes) or samples
#' @export
#' @docType methods
#' @rdname calcCoGAPSStat-methods
#'
#' @description calculates a statistic to determine if a pattern is enriched in a
#' a particular set of measurements or samples.
#' @param object an object of type CogapsResult
#' @param sets list of sets of measurements/samples
#' @param whichMatrix either "featureLoadings" or "sampleFactors" indicating which matrix
#' to calculate the statistics for
#' @param numPerm number of permutations to use when calculatin p-value
#' @param ... handles old arguments for backwards compatibility
#' @return gene set statistics for each column of A
setGeneric("calcCoGAPSStat", function(object, sets=NULL, whichMatrix='featureLoadings',
numPerm=1000, ...)
    {standardGeneric("calcCoGAPSStat")})

#' probability gene belongs in gene set
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
#' @param axis either 1 or 2, specifying if pattern markers should be calculated using
#' the rows of the data (1) or the columns of the data (2)
#' @return By default a non-overlapping list of genes associated with each
#' \code{lp}.
#' @examples
#' data(GIST)
#' pm <- patternMarkers(GIST.result)
setGeneric("patternMarkers", function(object, threshold="all", lp=NA, axis=1) standardGeneric("patternMarkers"))

#' MANOVA statistical test for patterns between sample groups
#' @export
#' @docType methods
#' @rdname MANOVA-methods
#' @description MANOVA statistical test--wraps base R manova
#' @param interestedVariables study design for manova
#' @param object CogapsResult object
#' @return list of manova fit results
setGeneric("MANOVA", function(interestedVariables, object) standardGeneric("MANOVA"))


#' save CoGAPS Result object as a set of csvs to directory 
#' see fromCSV
#' @export
#' @docType methods
#' @rdname toCSV-methods
#' @description save as csv
#' @param object CogapsResult object
#' @param save_location directory to write to
#' @return none
setGeneric("toCSV", function(object, save_location=".") standardGeneric("toCSV"))

#' read CoGAPS Result object from a directory with a set of csvs 
#' see toCSV
#' @export
#' @docType methods
#' @rdname fromCSV-methods
#' @description save as csv
#' @param save_location directory to read from
#' @return CogapsResult object
setGeneric("fromCSV", function(save_location=".") standardGeneric("fromCSV"))

