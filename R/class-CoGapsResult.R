setClass("CoGapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix"
))

setMethod("initialize", "CoGapsResult",
    function(.Object, in_Amean, in_Asd, in_Pmean, in_Psd)
    {
        # DONE?(Hyejune) store all input parameters into the
        # the object variables e.g. -- .Object@Amean <- in_Amean    
        .Object@Amean <- in_Amean
        .Object@Asd <- in_Asd
        .Object@Pmean <- in_Pmean
        .Object@Psd <- in_Psd
        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

# TODO(Hyejune, Tom) Override the plot function for CoGapsResult

setMethod('show', signature('CoGAPSResult'),
          function(object)
          {
            print(paste("Amean has", nrow(Amean), "rows and", ncol(Amean), "columns"))
            print(paste("Pmean has", nrow(Pmean), "rows and", ncol(Pmean), "columns"))
            
          }
)

setMethod('plot', signature('CoGAPSResult'),
          function(Pmean)
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
              
          }
)
