#' Plot the P Matrix
#'
#' @details plots the P matrix in a line plot with error bars
#' @param Pmean matrix of mean values of P
#' @param Psd matrix of standard deviation values of P
#' @export
plotP <- function(Pmean, Psd)
{
    Nfactor=dim(Pmean)[1]
    Nobs=dim(Pmean)[2]
    RowP <- 1:Nobs
    colors <- rainbow(Nfactor)
    ylimits <- c(0,(max(Pmean + Psd)*1.05))

    plotCI(x=RowP, y=Pmean[1,], col=colors[1],uiw=Psd[1,],
        ylim=ylimits,type='l',ylab="Relative Amplitude")

    for (i in 2:Nfactor)
    {
        points(RowP, Pmean[i,], col=colors[i], pch=i)
        plotCI(RowP, y=Pmean[i,], col=colors[i], uiw=Psd[i,],
            add=TRUE)
    }

    legend("bottom", paste("Pattern", 1:Nfactor, sep = ""),
        pch = 1:Nfactor, lty=1,cex=0.8, col=colors,bty="y",ncol=5)
}
