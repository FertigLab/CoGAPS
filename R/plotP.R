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
    RowP <- 1:Nobs
    colors <- rainbow(Nfactor)
    ylimits <- c(0,(max(PMean_Mat + P_SD)*1.05))

    plotCI(x=RowP, y=PMean_Mat[1,], col=colors[1],uiw=P_SD[1,],
          ylim=ylimits,type='l',ylab="Relative Amplitude")

    for (i in 2:Nfactor) {
        points(RowP, PMean_Mat[i,], col=colors[i], pch=i)
        plotCI(RowP, y=PMean_Mat[i,], col=colors[i], uiw=P_SD[i,],
            add=TRUE)
    }

    legend("bottom", paste("Pattern", 1:Nfactor, sep = ""),
         pch = 1:Nfactor, lty=1,cex=0.8,
        col=colors,bty="y",ncol=5)

}
