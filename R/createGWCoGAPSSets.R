#' Create Gene Sets for GWCoGAPS
#'
#' @details factors whole genome data into randomly generated sets for indexing
#'
#' @param D data matrix
#' @param S uncertainty matrix
#' @param nSets number of sets to partition the data into
#' @param simulationName name used to identify files created by this simulation
#' @return simulationName used to identify saved files
#' @examples
#' data(SimpSim)
#' createGWCoGAPSSets(SimpSim.D, SimpSim.S, nSets=2, "example")
#' @export
createGWCoGAPSSets <- function(D, S, nSets, simulationName)
{
    # check gene names
    if (length(unique(rownames(D))) != length(rownames(D)))
    {
        warning("Gene identifiers not unique!")
    }

    # partition data by sampling random sets of genes
    genes <- 1:nrow(D)
    setSize <- floor(length(genes) / nSets)
    for (set in 1:nSets)
    {
        # sample gene names
        sampleSize <- ifelse(set == nSets, length(genes), setSize)
        geneSet <- sample(genes, sampleSize, replace=FALSE)
        genes <- genes[!(genes %in% geneSet)]

        # partition data
        sampleD <- D[geneSet,]
        sampleS <- S[geneSet,]
        save(sampleD, sampleS, file=paste(simulationName, "_partition_", set,
            ".RData", sep=""));
    }
    return(simulationName)
}
