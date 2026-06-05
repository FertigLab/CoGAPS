#!/usr/bin/env Rscript

# this script is not used anywhere in the CoGAPS package, it has been written 
# by codex AI to debug reproducibility in CoGAPS results across different
# compiler versions and platforms.

suppressPackageStartupMessages(library(CoGAPS))

parseArgs <- function(args)
{
    defaults <- list(
        out = "debug_compiler_reproducibility_fingerprint.tsv",
        seed = 42,
        nIterations = 100,
        nPatterns = 7,
        outputFrequency = 10,
        repeats = 1
    )

    for (arg in args)
    {
        parts <- strsplit(sub("^--", "", arg), "=", fixed=TRUE)[[1]]
        if (length(parts) != 2)
            stop("arguments must use --name=value format: ", arg)
        name <- parts[1]
        value <- parts[2]
        if (!name %in% names(defaults))
            stop("unknown argument: ", name)
        if (name == "out")
            defaults[[name]] <- value
        else
            defaults[[name]] <- as.numeric(value)
    }

    defaults
}

hashObject <- function(x)
{
    path <- tempfile(fileext=".rds")
    on.exit(unlink(path), add=TRUE)
    saveRDS(x, path, version=2)
    unname(tools::md5sum(path))
}

summarizeMatrix <- function(x)
{
    c(
        hash=hashObject(x),
        nrow=nrow(x),
        ncol=ncol(x),
        sum=signif(sum(x), 16),
        mean=signif(mean(x), 16),
        min=signif(min(x), 16),
        max=signif(max(x), 16)
    )
}

resultFingerprint <- function(result)
{
    featureLoadings <- summarizeMatrix(result@featureLoadings)
    sampleFactors <- summarizeMatrix(result@sampleFactors)
    loadingStdDev <- summarizeMatrix(result@loadingStdDev)
    factorStdDev <- summarizeMatrix(result@factorStdDev)

    c(
        featureLoadings_hash=unname(featureLoadings["hash"]),
        featureLoadings_sum=unname(featureLoadings["sum"]),
        sampleFactors_hash=unname(sampleFactors["hash"]),
        sampleFactors_sum=unname(sampleFactors["sum"]),
        loadingStdDev_hash=unname(loadingStdDev["hash"]),
        loadingStdDev_sum=unname(loadingStdDev["sum"]),
        factorStdDev_hash=unname(factorStdDev["hash"]),
        factorStdDev_sum=unname(factorStdDev["sum"]),
        atomsA_hash=hashObject(result@metadata$atomsA),
        atomsP_hash=hashObject(result@metadata$atomsP),
        chisq_hash=hashObject(result@metadata$chisq),
        meanChiSq=signif(result@metadata$meanChiSq, 16)
    )
}

runMode <- function(data, opts, modeName, sparseOptimization,
asynchronousUpdates, nThreads, repeatIndex)
{
    effectiveAsynchronousUpdates <- asynchronousUpdates
    effectiveNThreads <- nThreads
    if (!CoGAPS::compiledWithOpenMPSupport())
    {
        effectiveAsynchronousUpdates <- FALSE
        effectiveNThreads <- 1
    }

    result <- NULL
    suppressWarnings(
        suppressMessages(
            capture.output(
                result <- CoGAPS(
                    data,
                    nPatterns=opts$nPatterns,
                    nIterations=opts$nIterations,
                    outputFrequency=opts$outputFrequency,
                    seed=opts$seed,
                    messages=FALSE,
                    sparseOptimization=sparseOptimization,
                    asynchronousUpdates=asynchronousUpdates,
                    nThreads=nThreads
                )
            )
        )
    )

    c(
        mode=modeName,
        repeatIndex=repeatIndex,
        seed=opts$seed,
        nIterations=opts$nIterations,
        nPatterns=opts$nPatterns,
        requestedSparseOptimization=sparseOptimization,
        requestedAsynchronousUpdates=asynchronousUpdates,
        requestedNThreads=nThreads,
        effectiveAsynchronousUpdates=effectiveAsynchronousUpdates,
        effectiveNThreads=effectiveNThreads,
        resultFingerprint(result)
    )
}

opts <- parseArgs(commandArgs(trailingOnly=TRUE))

data(GIST)
data <- GIST.matrix

buildReportText <- gsub("[\r\n\t]+", " | ", CoGAPS::buildReport())
metadata <- c(
    packageVersion=as.character(utils::packageVersion("CoGAPS")),
    RVersion=R.version.string,
    platform=R.version$platform,
    openMP=CoGAPS::compiledWithOpenMPSupport(),
    buildReport=buildReportText,
    os=paste(Sys.info()[c("sysname", "release", "machine")], collapse=" ")
)

modes <- list(
    list(name="dense_sequential", sparseOptimization=FALSE,
        asynchronousUpdates=FALSE, nThreads=1),
    list(name="sparse_sequential", sparseOptimization=TRUE,
        asynchronousUpdates=FALSE, nThreads=1),
    list(name="dense_async_one_thread", sparseOptimization=FALSE,
        asynchronousUpdates=TRUE, nThreads=1),
    list(name="sparse_async_one_thread", sparseOptimization=TRUE,
        asynchronousUpdates=TRUE, nThreads=1),
    list(name="dense_async_two_threads", sparseOptimization=FALSE,
        asynchronousUpdates=TRUE, nThreads=2),
    list(name="sparse_async_two_threads", sparseOptimization=TRUE,
        asynchronousUpdates=TRUE, nThreads=2)
)

rows <- list()
for (repeatIndex in seq_len(opts$repeats))
{
    for (mode in modes)
    {
        fingerprint <- runMode(
            data=data,
            opts=opts,
            modeName=mode$name,
            sparseOptimization=mode$sparseOptimization,
            asynchronousUpdates=mode$asynchronousUpdates,
            nThreads=mode$nThreads,
            repeatIndex=repeatIndex
        )
        rows[[length(rows) + 1]] <- c(metadata, fingerprint)
    }
}

out <- as.data.frame(do.call(rbind, rows), stringsAsFactors=FALSE)
write.table(out, file=opts$out, quote=TRUE, sep="\t", row.names=FALSE)
cat("wrote compiler reproducibility fingerprint to", opts$out, "\n")
