#' get specified number of retina subsets
#' @export
#'
#' @description combines retina subsets from extdata directory
#' @param n number of subsets to use
#' @return matrix of RNA counts
#' @examples
#' retSubset <- getRetinaSubset()
#' dim(retSubset)
#' @importFrom rhdf5 h5read
getRetinaSubset <- function(n=1)
{
    if (!(n %in% 1:4))
        stop("invalid number of subsets requested")

    subset_1_path <- system.file("extdata/retina_subset_1.h5", package="CoGAPS")
    subset_2_path <- system.file("extdata/retina_subset_2.h5", package="CoGAPS")
    subset_3_path <- system.file("extdata/retina_subset_3.h5", package="CoGAPS")
    subset_4_path <- system.file("extdata/retina_subset_4.h5", package="CoGAPS")

    data <- rhdf5::h5read(subset_1_path, "counts")
    cNames <- rhdf5::h5read(subset_1_path, "cellNames")
    if (n > 1)
    {
        data <- cbind(data, rhdf5::h5read(subset_2_path, "counts"))
        cNames <- c(cNames, rhdf5::h5read(subset_2_path, "cellNames"))
    }
    if (n > 2)
    {
        data <- cbind(data, rhdf5::h5read(subset_3_path, "counts"))
        cNames <- c(cNames, rhdf5::h5read(subset_3_path, "cellNames"))
    }
    if (n > 3)
    {
        data <- cbind(data, rhdf5::h5read(subset_4_path, "counts"))
        cNames <- c(cNames, rhdf5::h5read(subset_4_path, "cellNames"))
    }

    colnames(data) <- cNames
    rownames(data) <- rhdf5::h5read(subset_1_path, "geneNames")
    return(data)    
}

#' wrapper around cat
#' @keywords internal
#'
#' @description cleans up message printing
#' @param allParams all cogaps parameters
#' @param ... arguments forwarded to cat
#' @return conditionally print message
gapsCat <- function(allParams, ...)
{
    if (allParams$messages)
        cat(...)
}

#' checks if file is supported
#' @keywords internal
#'
#' @param file path to file
#' @return TRUE if file is supported, FALSE if not
#' @importFrom tools file_ext
supported <- function(file)
{
    if (!is(file, "character"))
        return(FALSE)
    return(tools::file_ext(file) %in% c("tsv", "csv", "mtx", "gct"))
}

#' checks if file is rds format
#' @keywords internal
#'
#' @param file path to file
#' @return TRUE if file is .rds, FALSE if not
#' @importFrom tools file_ext
isRdsFile <- function(file)
{
    if (is.null(file))
        return(FALSE)
    if (length(file) == 0)
        return(FALSE)
    if (!is(file, "character"))
        return(FALSE)
    return(tools::file_ext(file) == "rds")
}

#' get input that might be an RDS file
#' @keywords internal
#'
#' @param input some user input
#' @return if input is an RDS file, read it - otherwise return input
getValueOrRds <- function(input)
{
    if (isRdsFile(input))
        return(readRDS(input))
    return(input)
}

#' get number of rows from supported file name or matrix
#' @keywords internal
#'
#' @param data either a file name or a matrix
#' @return number of rows
#' @importFrom data.table fread
#' @importFrom tools file_ext
nrowHelper <- function(data)
{
    if (is(data, "character"))
    {
        return(getFileInfo_cpp(data)[["dimensions"]][1])
    }
    return(nrow(data))
}

#' get number of columns from supported file name or matrix
#' @keywords internal
#'
#' @param data either a file name or a matrix
#' @return number of columns
#' @importFrom data.table fread
#' @importFrom tools file_ext
ncolHelper <- function(data)
{
    if (is(data, "character"))
    {
        return(getFileInfo_cpp(data)[["dimensions"]][2])
    }
    return(ncol(data))
}

#' extract gene names from data
#' @keywords internal
#' @return vector of gene names
getGeneNames <- function(data, transpose)
{
    if (transpose)
        return(getSampleNames(data, FALSE))
    if (is(data, "character"))
        names <- getFileInfo_cpp(data)[["rowNames"]]
    else
        names <- rownames(data)
    if (is.null(names) | length(names) == 0)
        return(paste("Gene", 1:nrowHelper(data), sep="_"))
    return(names)
}

#' extract sample names from data
#' @keywords internal
#' @return vector of sample names
getSampleNames <- function(data, transpose)
{
    if (transpose)
        return(getGeneNames(data, FALSE))
    if (is(data, "character"))
        names <- getFileInfo_cpp(data)[["colNames"]]
    else
        names <- colnames(data)
    if (is.null(names) | length(names) == 0)
        return(paste("Sample", 1:ncolHelper(data), sep="_"))
    return(names)
}

#' write start up message
#' @keywords internal
#'
#' @param data data set
#' @param allParams list of all parameters
#' @return message displayed to screen
#' @importFrom methods show
startupMessage <- function(data, allParams)
{
    nGenes <- ifelse(allParams$transposeData, ncolHelper(data), nrowHelper(data))
    nSamples <- ifelse(allParams$transposeData, nrowHelper(data), ncolHelper(data))

    dist_message <- "Standard"
    if (!is.null(allParams$gaps@distributed))
        dist_message <- allParams$gaps@distributed

    cat("\nThis is CoGAPS version", as.character(packageVersion("CoGAPS")), "\n")
    cat("Running", dist_message, "CoGAPS on", allParams$dataName,
        paste("(", nGenes, " genes and ", nSamples, " samples)", sep=""))

    if (allParams$messages)
    {
        cat(" with parameters:\n\n")
        methods::show(allParams$gaps)
    }
    cat("\n")
}

#' parse parameters passed through the ... variable
#' @keywords internal
#'
#' @param allParams list of all parameters
#' @param extraParams list of parameters in ...
#' @return allParams with any valid parameters in extraParams added
#' @note will halt with an error if any parameters in extraParams are invalid
#' @importFrom methods slotNames
parseExtraParams <- function(allParams, extraParams)
{
    # parse direct params
    for (s in slotNames(allParams$gaps))
    {
        if (!is.null(extraParams[[s]]))
        {
            allParams$gaps <- setParam(allParams$gaps, s, extraParams[[s]])
            extraParams[[s]] <- NULL
        }
    }

    # check for unrecognized options
    if (length(extraParams) > 0)
        stop(paste("unrecognized argument:", names(extraParams)[1]))

    return(allParams)
}

## TODO these checks should be in the C++ code so that file names are checked
## just as much as R variables
#' check that provided data is valid
#' @keywords internal
#'
#' @param data data matrix
#' @param uncertainty uncertainty matrix, can be null
#' @param params CogapsParams object
#' @return throws an error if data has problems
checkDataMatrix <- function(data, uncertainty, params)
{
    if (any(is.na(data)))
        stop("NA values in data")
    if (!all(apply(data, 2, is.numeric)))
        stop("data is not numeric")
    if (sum(data < 0) > 0 | sum(uncertainty < 0) > 0)
        stop("negative values in data and/or uncertainty matrix")
    if (nrow(data) <= params@nPatterns | ncol(data) <= params@nPatterns)
        stop("nPatterns must be less than dimensions of data")
    if (sum(uncertainty < 1e-5) > 0)
        warning("small values in uncertainty matrix detected")
}

#' check that all inputs are valid
#' @keywords internal
#'
#' @param data data matrix
#' @param uncertainty uncertainty matrix, can be null
#' @param allParams list of all parameters
#' @return throws an error if inputs are invalid
checkInputs <- function(data, uncertainty, allParams)
{
    if (is(data, "character") & !is.null(uncertainty) & !is(uncertainty, "character"))
        stop("uncertainty must be same data type as data (file name)")
    if (is(uncertainty, "character") & !supported(uncertainty))
        stop("unsupported file extension for uncertainty")
    if (!is(data, "character") & !is.null(uncertainty) & !is(uncertainty, "matrix"))
        stop("uncertainty must be a matrix unless data is a file path")
    if (!is.null(uncertainty) & allParams$gaps@sparseOptimization)
        stop("must use default uncertainty when enabling sparseOptimization")
    if (!is.null(allParams$checkpointInFile) & !CoGAPS::checkpointsEnabled())
        stop("CoGAPS was built with checkpoints disabled")

    if (!is.null(allParams$gaps@distributed))
    {
        if (allParams$gaps@distributed == "single-cell" & !allParams$gaps@singleCell)
            warning("running single-cell CoGAPS with singleCell=FALSE")
        if (allParams$nThreads > 1)
            stop("can't run multi-threaded and distributed CoGAPS at the same time")
        if (!is.null(allParams$checkpointInFile))
            stop("checkpoints not supported for distributed cogaps")
        if (!is(data, "character"))
            warning("running distributed cogaps without mtx/tsv/csv/gct data")
    }

    if (!is(data, "character"))
        checkDataMatrix(data, uncertainty, allParams$gaps)
}

#' extracts gene/sample names from the data
#' @keywords internal
#'
#' @param data data matrix
#' @param allParams list of all parameters
#' @return list of all parameters with added gene names
getDimNames <- function(data, allParams)
{
    # get user supplied names
    geneNames <- allParams$gaps@geneNames
    sampleNames <- allParams$gaps@sampleNames

    # if user didn't supply any names, pull from data set or use default labels
    if (is.null(allParams$gaps@geneNames))
        geneNames <- getGeneNames(data, allParams$transposeData)
    if (is.null(allParams$gaps@sampleNames))
        sampleNames <- getSampleNames(data, allParams$transposeData)

    # get the number of genes/samples
    nGenes <- ifelse(allParams$transposeData, ncolHelper(data), nrowHelper(data))
    nSamples <- ifelse(allParams$transposeData, nrowHelper(data), ncolHelper(data))

    # handle any subsetting
    if (allParams$gaps@subsetDim == 1)
    {
        nGenes <- length(allParams$gaps@subsetIndices)
        geneNames <- geneNames[allParams$gaps@subsetIndices]
    }
    else if (allParams$gaps@subsetDim == 2)
    {
        nSamples <- length(allParams$gaps@subsetIndices)
        sampleNames <- sampleNames[allParams$gaps@subsetIndices]
    }    

    # check that names align with expected number of genes/samples
    if (length(geneNames) != nGenes)
        stop(length(geneNames), " != ", nGenes, " incorrect number of gene names given")
    if (length(sampleNames) != nSamples)
        stop(length(sampleNames), " != ", nSamples, " incorrect number of sample names given")

    # store processed gene/sample names directly in allParams list
    # this is an important distinction - allParams@gaps contains the
    # gene/sample names originally passed by the user, allParams contains
    # the procseed gene/sample names to be used when labeling the result
    allParams$geneNames <- geneNames
    allParams$sampleNames <- sampleNames
    return(allParams)
}

#' convert any acceptable data input to a numeric matrix
#' @keywords internal
#'
#' @description convert supported R objects containing the data to a
#' numeric matrix, if data is a file name do nothing. Exits with an error
#' if data is not a supported type.
#' @param data data input
#' @return data matrix
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
convertDataToMatrix <- function(data)
{
    if (is(data, "character") & !supported(data))
        stop("unsupported file extension for data")
    else if (is(data, "matrix") | is(data, "character"))
        return(data)
    else if (is(data, "data.frame"))
        return(data.matrix(data))
    else if (is(data, "SummarizedExperiment"))
        return(SummarizedExperiment::assay(data, "counts"))
    else if (is(data, "SingleCellExperiment"))
        return(SummarizedExperiment::assay(data, "counts"))
    else
        stop("unsupported data type")
}