#' get specified number of retina subsets
#' @export
#'
#' @description combines retina subsets from extdata directory
#' @param n number of subsets to use
#' @return matrix of RNA counts
#' @importFrom rhdf5 h5read
getRetinaSubset <- function(n=1)
{
    if (!(n %in% 1:4))
    {
        stop("invalid number of subsets requested")
    }

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
        return(switch(tools::file_ext(data),
            "csv" = nrow(data.table::fread(data, select=1)),
            "tsv" = nrow(data.table::fread(data, select=1)),
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,1]),
            "gct" = as.numeric(strsplit(as.matrix(data.table::fread(data, nrows=1, sep='\t')), "\\s+")[[1]][1])
        ))
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
        return(switch(tools::file_ext(data),
            "csv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "tsv" = ncol(data.table::fread(data, nrows=1)) - 1,
            "mtx" = as.numeric(data.table::fread(data, nrows=1, fill=TRUE)[1,2]),
            "gct" = as.numeric(strsplit(as.matrix(data.table::fread(data, nrows=1, sep='\t')), "\\s+")[[1]][2])
        ))
    }
    return(ncol(data))
}

#' extract gene names from data
#' @keywords internal
#' @return vector of gene names
getGeneNames <- function(data, transpose)
{
    if (transpose)
    {
        return(getSampleNames(data, FALSE))
    }

    names <- NULL
    if (is(data, "character"))
    {
        names <- switch(tools::file_ext(data),
            "csv" = as.matrix(data.table::fread(data, select=1))[,1],
            "tsv" = as.matrix(data.table::fread(data, select=1))[,1],
            "gct" = suppressWarnings(gsub("\"", "", as.matrix(data.table::fread(data, select=1))))
        )
    }
    else if (is(data, "matrix") | is(data, "data.frame"))
    {
        names <- rownames(data)
    }

    if (is.null(names))
    {
        nGenes <- nrowHelper(data)
        return(paste("Gene", 1:nGenes, sep="_"))
    }
    return(names)
}

#' extract sample names from data
#' @keywords internal
#' @return vector of sample names
getSampleNames <- function(data, transpose)
{
    if (transpose)
    {
        return(getGeneNames(data, FALSE))
    }

    names <- NULL
    if (is(data, "character"))
    {
        names <- switch(tools::file_ext(data),
            "csv" = colnames(data.table::fread(data, nrows=1))[-1],
            "tsv" = colnames(data.table::fread(data, nrows=1))[-1],
            "gct" = suppressWarnings(colnames(data.table::fread(data, skip=2, nrows=1))[-1:-2])
        )
    }
    else if (is(data, "matrix") | is(data, "data.frame"))
    {
        names <- colnames(data)
    }

    if (is.null(names))
    {
        nSamples <- ncolHelper(data)
        return(paste("Sample", 1:nSamples, sep="_"))
    }
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
    cat("Running", dist_message, "CoGAPS on", nGenes, "genes and",
        nSamples, "samples")

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

#' check that provided data is valid
#' @keywords internal
#'
#' @param data data matrix
#' @param uncertainty uncertainty matrix, can be null
#' @param params CogapsParams object
#' @return throws an error if data has problems
checkDataMatrix <- function(data, uncertainty, params)
{
    if (sum(data < 0) > 0 | sum(uncertainty < 0) > 0)
        stop("negative values in data and/or uncertainty matrix")
    if (nrow(data) <= params@nPatterns | ncol(data) <= params@nPatterns)
        stop("nPatterns must be less than dimensions of data")
    if (sum(uncertainty < 1e-5) > 0)
        warning("small values in uncertainty matrix detected")
}