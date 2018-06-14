setClass("CoGapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix"
))

setMethod("initialize", "CoGapsResult",
    function(.Object, in_Amean, in_Asd, in_Pmean, in_Psd)
    {
        # TODO(Hyejune) store all input parameters into the
        # the object variables e.g. -- .Object@Amean <- in_Amean    

        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)

# TODO(Hyejune, Tom) Override the plot function for CoGapsResult

