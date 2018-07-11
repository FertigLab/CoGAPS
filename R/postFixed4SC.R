#' Post Processing of Parallel Output
#'
#' @param AP.fixed output of parallel gapsMapRun calls with same FP
#' @param setAs data.frame with rows giving fixed patterns for A used as input
#' for gapsMapRun
#' @return list of two data.frames containing the A matrix estimates or their
#' corresponding standard deviations from output of parallel CoGAPS
postFixed4SC <- function(AP.fixed, setAs)
{
    ASummary <- AP.fixed[[1]]$Amean
    PSummary <- do.call(cbind,lapply(AP.fixed, function(x) t(x@sampleFactors)))
    Psd <- do.call(cbind,lapply(AP.fixed, function(x) t(x@sampleStdDev)))

    Amax <- apply(ASummary,2,max)
    Aneu <- sweep(ASummary,2,Amax,FUN="/")
    Pneu <- sweep(PSummary,1,Amax,FUN="*")

    X <- apply(Aneu,2,range)
    Y <- apply(setAs,2,range)
    colnames(X) <- colnames(Y)
    if (all.equal(X,Y,tolerance=0.01) != TRUE)
        warning("Patterns do not match fixed values.")

    Ps4fixAs<-list("P"=Pneu,"Psd"=Psd)
    return(Ps4fixAs)
}
