#' CogapsResult
#' @export
#'
#' @description Contains all output from Cogaps run
setClass("CogapsResult", slots=c(
    Amean = "matrix",
    Asd = "matrix",
    Pmean = "matrix",
    Psd = "matrix",
    seed = "numeric",
    meanChiSq = "numeric",
    diagnostics = "list"
))

#' Constructor for CogapsResult
#' @return initialized CogapsResult object
#' @importFrom methods callNextMethod
setMethod("initialize", "CogapsResult",
    function(.Object, ...)
    {
        .Object <- callNextMethod(.Object, ...)
        .Object
    }
)
