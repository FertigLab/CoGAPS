# plotP: function to plot Rows of P with error bars
# History: v1.0 - JF 7/2014

# Inputs: PMean_Mat - matrix with rows to plot
#         P_SD - uncertainties on PMean_Mat for error bars

# Output: screen plot of rows of PMean_Mat

#'\code{plotP} plots the P matrix in a line plot
#' with error bars
#'
#'@param PMean_Mat matrix of mean values of P
#'@param P_SD matrix of standard deviation values of P
#'@export


plotP<-function(PMean_Mat, P_SD)  {


    Nfactor=dim(PMean_Mat)[1]
    Nobs=dim(PMean_Mat)[2]
    RowP <- rep(1:Nobs,times=Nfactor)
    max <- nrow(PMean_Mat)

    layout(rbind(1,2), heights=c(9,1)) 
    matplot(t(PMean_Mat), type="b", pch = 1:max,cex=0.6, lty=rep(1,max), col = rainbow(max),
        main="Patterns from P Matrix",ylab="Amplitude",
        cex.axis=0.8,cex.main=1.4,cex.lab=1.3) 
    plotCI(x=RowP, y= as.vector(t(PMean_Mat)), 
        uiw=c(t(matrix(P_SD,nrow=Nfactor,ncol=Nobs,byrow=FALSE))), 
        add=T,col = rep(rainbow(ncol(t(PMean_Mat))),
        each=nrow(t(PMean_Mat)))) 
    par(mar=c(0, 0, 0, 0))
    plot.new()
    legend("top", paste("Pattern", 1:max, sep = ""), pch = 1:max, lty=1,cex=0.8,
        col=rainbow(max),bty="y",ncol=5)

}
