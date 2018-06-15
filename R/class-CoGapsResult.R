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

