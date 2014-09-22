#residuals: function to create a heatmap for the residuals
# History: v1.0 - JS, KC 7/2014

# Inputs: AMean_Mat - mean of A matrix
#         PMean_Mat - mean of P matrix
#         D - input data matrix
#         S - input standard deviation matrix

# Output: screen plot of residuals

#'\code{residuals} calculate residuals and produce heatmap
#'
#'@param AMean_Mat matrix of mean values for A from GAPS
#'@param PMean_Mat matrix of mean values for P from GAPS
#'@param D original data matrix run through GAPS
#'@param S original standard deviation matrix run through GAPS
#'@export


residuals=function(AMean_Mat, PMean_Mat, D, S)    
{

    M_Mean <- AMean_Mat%*%PMean_Mat
    Resid_M_Mean<-as.matrix((D - M_Mean)/S)
    colnames(Resid_M_Mean) <- colnames(D)
    rownames(Resid_M_Mean) <- rownames(D)

    scaledRdYlBu <- colorRampPalette(brewer.pal(9,"RdYlBu"))(100)
    heatmap.2(Resid_M_Mean, Rowv = FALSE, Colv = FALSE,dendrogram="none",
        scale="none",col = scaledRdYlBu, trace="none",density.info="none",
        cexCol=1.33,srtCol=45,lmat=rbind(c(0, 3),c(2,1),c(0,4) ),
        lwid=c(1,10),lhei=c(1, 4, 1.2 ), main="Heatmap of Residuals")
}
