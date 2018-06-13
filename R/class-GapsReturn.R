GapsReturn <- setClass("GapsReturn", slots=c(
  Amean = "matrix",
  Asd = "matrix",
  Pmean = "matrix",
  Psd = "matrix"
))

setMethod("initialize", "GapsReturn",
          function(returnList, returnObject, ...)
          {
            returnList <- cogaps_cpp(...)
            returnObject <- new('GapsReturn')
            returnObject@Amean <- returnList$Amean
            returnObject@Asd <- returnList$Asd
            returnObject@Pmean <- returnList$Pmean
            returnObject@Psd <- returnList$Psd
            return(returnObject)
          }
)

