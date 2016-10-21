#' plotPatternMarkers
#'
#' @param data the dataset from which the patterns where generated
#' @param patternMarkers the list of genes generated from the patternMarkers function
#' @param patternPalette a vector indicating what color should be used for each pattern
#' @param sampleNames names of the samples to use for labeling 
#' @param samplePalette  a vector indicating what color should be used for each sample
#' @param colDenogram logical indicating whether to display sample denogram
#' @param heatmapCol pallelet giving color scheme for heatmap
#' @param scale character indicating if the values should be centered and scaled in either 
#' the row direction or the column direction, or none. The default is "row". 
#' @param ... additional graphical parameters to be passed to \code{heatmap.2}
#' @export
#' @return heatmap of the \code{data} values for the \code{patternMarkers}
#' @seealso  \code{\link{heatmap.2}}
#' @examples \dontrun{
#' plotPatternMarkers(data=p,patternMarkers=PatternMarkers,patternPalette=NA,sampleNames=pd$sample,
#' samplePalette=pd$color,colDenogram=TRUE,heatmapCol="bluered", scale='row')
#'}

plotPatternMarkers <- function(
    data=NA,# the dataset from which the patterns where generated
    patternMarkers=NA,# the list of genes generated from the patternMarkers function
    patternPalette=NA,# a vector indicating what color should be used for each pattern
    sampleNames=NA,#names of the samples to use for labeling
    samplePalette=NULL,# a vector indicating what color should be used for each sample
    colDenogram=TRUE,# logical indicating whether to
    heatmapCol="bluered",
    scale='row',
    ...){

if(is.null(samplePalette)){samplePalette<-rep("black",dim(data)[2])}

### coloring of genes
patternCols <- rep(NA, length(unlist(patternMarkers)))
names(patternCols) <- unlist(patternMarkers)
for (i in 1:length(patternMarkers)) {
    patternCols[patternMarkers[[i]]] <- patternPalette[i]
}

heatmap.2(as.matrix(data[unlist(patternMarkers),],...),
          symkey=F,  symm=F, 
          scale=scale,
          col=heatmapCol,
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          hclustfun=function(x) hclust(x,method="average"),
          Rowv=F,Colv=colDenogram,trace='none',
          RowSideColors = as.character(patternCols[unlist(patternMarkers)]),
          labCol= sampleNames,
          cexCol=.8,
          ColSideColors = as.character(samplePalette),
          rowsep = cumsum(sapply(patternMarkers,length)))
}
