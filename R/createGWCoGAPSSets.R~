#' createGWCoGAPSSets
#'
#'\code{createGWCoGAPSSets} factors whole genome data into randomly generated sets for indexing;
#'
#'@param D data matrix with unique rownames
#'@param nSets number of sets for parallelization
#'@param outRDA name of output file
#'@param keep logical indicating whether or not to save gene set list. Default is TRUE.
#'@export
#'@return list with randomly generated sets of genes from whole genome data
#'@examples \dontrun{
#'createGWCoGAPSSet(D,nSets=nSets)
#'}
#'

createGWCoGAPSSets<-function(data=D, #data matrix with unique rownames
	nSets=nSets, #number of sets for parallelization
	outRDA="GenesInCoGAPSSets.Rda", #name of output file
	keep=TRUE #logical indicating whether or not to save gene set list. Default is TRUE.
	){
genes=rownames(data)
setSize=floor(length(genes)/nSets)
genesInSets <- list()
for (set in 1:nSets) {
  if(set!=nSets){genesInSets[[set]] <- sample(genes,setSize)}
  if(set==nSets){genesInSets[[set]] <- genes}
  genes=genes[!genes%in%genesInSets[[set]]]
}
if(!identical(sort(unlist(genesInSets)),sort(rownames(data)))){print("Gene identifiers not unique!")}
if(keep==TRUE){save(list=c('genesInSets'),file=outRDA)}
return(genesInSets)
}


