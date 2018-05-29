#' Plot the P Matrix
#'
#' @details plots the P matrix in a line plot with error bars
#' @param Pmean matrix of mean values of P
#' @param Psd matrix of standard deviation values of P
#' @return plot
#' @examples
#' data(SimpSim)
#' plotP(SimpSim.result$Pmean, SimpSim.result$Psd)
#' @importFrom gplots plotCI
#' @export
plotP <- function(Pmean, Psd=NULL)
{
    colors <- rainbow(nrow(Pmean))
    xlimits <- c(0, ncol(Pmean) + 1)
    ylimits <- c(0, (max(Pmean) * 1.05))

    plot(NULL, xlim=xlimits, ylim=ylimits, ylab="Relative Amplitude")

    for (i in 1:nrow(Pmean))
    {
        lines(x=1:ncol(Pmean), y=Pmean[i,], col=colors[i])
        points(x=1:ncol(Pmean), y=Pmean[i,], col=colors[i], pch=i)
    }
    
    legend("bottom", paste("Pattern", 1:nrow(Pmean), sep = ""),
        pch = 1:nrow(Pmean), lty=1, cex=0.8, col=colors, bty="y", ncol=5)

    #Nfactor <- nrow(Pmean)
    #Nobs <- ncol(Pmean)
    #RowP <- 1:Nobs
    #colors <- rainbow(Nfactor)
    #ylimits <- c(0,(max(Pmean + Psd)*1.05))
    #
    #plotCI(x=RowP, y=Pmean[1,], col=colors[1], uiw=Psd[1,],
    #    ylim=ylimits, type='l', ylab="Relative Amplitude")
    #
    #for (i in 2:Nfactor)
    #{
    #    #points(RowP, Pmean[i,], col=colors[i], pch=i)
    #    plotCI(RowP, y=Pmean[i,], col=colors[i], uiw=Psd[i,],
    #        add=TRUE)
    #}
    #
    #legend("bottom", paste("Pattern", 1:Nfactor, sep = ""),
    #    pch = 1:Nfactor, lty=1,cex=0.8, col=colors,bty="y",ncol=5)
}
