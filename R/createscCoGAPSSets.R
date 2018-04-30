#' Create Gene Sets for scCoGAPS
#' @export
#'
#' @description factors whole genome data into randomly generated sets for indexing
#' @param D data matrix
#' @param nSets number of sets to partition the data into
#' @param simulationName name used to identify files created by this simulation
#' @param samplingRatio vector of relative quantities to use for sampling celltypes
#' @param anotionObj vector of same length as number of columns of D 
#' @param path character string indicating were to save resulting data objects. default is current working dir
#' @return simulationName used to identify saved files
#' @examples
#' data(SimpSim)
#' createscCoGAPSSets(SimpSim.D, nSets=2, simulationName="example")
createscCoGAPSSets <- function(D, nSets, simulationName, samplingRatio=NULL,
path="", anotionObj=NULL)
{
    # check gene names
    if (length(unique(colnames(D))) != length(colnames(D)))
    {
        warning("Cell identifiers not unique!")
    }

    # partition data by sampling random sets of cells
    cells <- 1:ncol(D)
    setSize <- floor(length(cells) / nSets)
    for (set in 1:nSets)
    {
        if (is.null(samplingRatio))
        {
            # sample cell names
            sampleSize <- ifelse(set == nSets, length(cells), setSize)
            cellset <- sample(cells, sampleSize, replace=FALSE)
            cells <- cells[!(cells %in% cellset)]
        }
        else
        {
            if (length(unique(anotionObj)) != length(samplingRatio))
            {
                warning("Not all celltypes will be sampled from.")
            }
            ct.indx <- lapply(unique(anotionObj), function(x) which(anotionObj == x))
            cellset <- lapply(unique(anotionObj), function(x)
                sample(colnames(D)[ct.indx[[x]]], samplingRatio[x],replace=TRUE))
        }

        # partition data
        sampleD <- D[,cellset]
        #log transform 
        sampleD <- log2(sampleD+1)
        # generate S
        sampleS <- pmax(.1*sampleD, .1)
        save(sampleD, sampleS, file=paste0(path,simulationName, "_partition_",
            set,".RData"));
    }
    return(simulationName)
}
