#' Post Processing of Parallel Output
#'
#' @param AP.fixed output of parallel gapsMapRun calls with same FP
#' @param setPs data.frame with rows giving fixed patterns for P used as input
#' for gapsMapRun
#' @return list of two data.frames containing the A matrix estimates or their
#' corresponding standard deviations from output of parallel CoGAPS
postFixed4SC <- function(AP.fixed, setAs)
{
    ASummary <- AP.fixed[[1]]$Amean


    PSummary <- do.call(cbind,lapply(AP.fixed, function(x) x$Pmean))
    Psd <- do.call(cbind,lapply(AP.fixed, function(x) x$Psd))

    Pmax <- apply(PSummary,2,max)
    Pneu <- sweep(PSummary,2,Pmax,FUN="/")
    Aneu <- sweep(ASummary,1,Pmax,FUN="*")

    X <- apply(Pneu,2,range)
    Y <- apply(setPs,2,range)
    colnames(X) <- colnames(Y)
    if (all.equal(X,Y,tolerance=0.01) != TRUE)
        warning("Patterns do not match fixed values.")

    Ps4fixAs<-list("P"=Pneu,"Psd"=Psd)
    return(Ps4fixAs)
}
