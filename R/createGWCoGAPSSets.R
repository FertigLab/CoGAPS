#' createGWCoGAPSSets
#'
#'\code{createGWCoGAPSSets} factors whole genome data into randomly generated sets for indexing;
#'
#'@param data data matrix with unique rownames
#'@param nSets number of sets for parallelization
#'@param outRDA name of output file
#'@param keep logical indicating whether or not to save gene set list. Default is TRUE.
#'@param path character string indicating were to save resulting data objects
#'@export
#'@return list with randomly generated sets of genes from whole genome data
#'@examples \dontrun{
#'createGWCoGAPSSet(D,nSets=nSets)
#'}
#'

createGWCoGAPSSets <- function(D, S, nSets, simulationName, path="")
{
    # check gene names
    if (length(unique(colnames(D))) != length(colnames(D)))
    {
        warning("Cell identifiers not unique!")
    }

    # partition data by sampling random sets of cells
    genes <- 1:nrow(D)
    setSize <- floor(length(genes) / nSets)
    for (set in 1:nSets)
    {
        
        # sample genes
            sampleSize <- ifelse(set == nSets, length(genes), setSize)
            geneset <- sample(genes, sampleSize, replace=FALSE)
            genes <- genes[!(genes %in% geneset)]
        # partition data
        sampleD <- D[geneset,]
        sampleS <- S[geneset,]
        save(sampleD, sampleS, file=paste(path, simulationName, "_partition_", set,
            ".RData", sep=""));
    }
    return(simulationName)
}
