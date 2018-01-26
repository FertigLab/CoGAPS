#' Binary Heatmap for Standardized A Matrix
#'
#' @details creates a binarized heatmap of the A matrix
#'  in which the value is 1 if the value in Amean is greater than
#'  threshold * Asd and 0 otherwise
#' @param Amean the mean estimate for the A matrix
#' @param Asd the standard deviations on Amean
#' @param threshold the number of standard deviations above zero
#'  that an element of Amean must be to get a value of 1
#' @examples
#' # Load the outputs from gapsRun
#' data('results')
#' # Run binaryA with the correct arguments from 'results'
#' binaryA(results$Amean,results$Asd,threshold=3)
#' @export
binaryA <-function(Amean, Asd, threshold=3)
{
    BinA_Map <- ifelse(Amean/Asd > threshold,1,0)
    colnames(BinA_Map) <- colnames(BinA_Map, do.NULL=FALSE, prefix = "Pattern ")
    rownames(BinA_Map) <- rep(" ",nrow(BinA_Map))

    heatmap.2(BinA_Map, Rowv = FALSE, Colv = FALSE,dendrogram="none",
        scale="none",col = brewer.pal(3,"Blues"), trace="none",
        density.info="none",cexCol=1.3,srtCol=45,
        lmat=rbind(c(0, 3),c(2,1),c(0,4) ),
        lwid=c(1,10),lhei=c(1, 4, 1.2 ),
        main="Heatmap of Standardized A Matrix")
    mtext(paste("(Threshold = ", threshold, ")"))
}
