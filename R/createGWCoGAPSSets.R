#' Create Gene Sets for GWCoGAPS
#'
#' @details factors whole genome data into randomly generated sets for indexing
#'
#' @param data data matrix with unique rownames
#' @param nSets number of sets for parallelization
#' @param outRDA name of output file
#' @param keep logical indicating whether or not to save gene set list.
#' @return list with randomly generated sets of genes from whole genome data
#' @examples \dontrun{
#'createGWCoGAPSSet(D,nSets=nSets)
#'}
#' @export
createGWCoGAPSSets<-function(data=D, nSets=nSets,
outRDA="GenesInCoGAPSSets.Rda", keep=TRUE)
{
    genes=rownames(data)
    setSize=floor(length(genes)/nSets)
    genesInSets <- list()
    for (set in 1:nSets)
    {
        if(set!=nSets)
            genesInSets[[set]] <- sample(genes,setSize)
        if(set==nSets)
            genesInSets[[set]] <- genes
        genes=genes[!genes%in%genesInSets[[set]]]
    }
    if (!identical(sort(unlist(genesInSets)),sort(rownames(data))))
        message("Gene identifiers not unique!")
    if (keep==TRUE)
        save(list=c('genesInSets'),file=outRDA)
    return(genesInSets)
}

