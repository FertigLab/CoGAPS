#' Post Processing of Parallel Output
#'
#' @param AP.fixed output of parallel gapsMapRun calls with same FP
#' @param setValues data.frame with rows giving fixed patterns for P used as input
#' for gapsMapRun
#' @param setMatrix which matrix, A or P
#' @return list of two data.frames containing the A matrix estimates or their
#' corresponding standard deviations from output of parallel CoGAPS
postFixed4Parallel <- function(AP.fixed, setValues, setMatrix="P")
{
    if (setMatrix=="P")
    {
        ASummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Amean))
        Asd <- do.call(rbind,lapply(AP.fixed, function(x) x$Asd))
        #PSummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Pmean))
        PSummary <- AP.fixed[[1]]$Pmean

        Pmax <- apply(PSummary,1,max)
        Pneu <- sweep(PSummary,1,Pmax,FUN="/")
        Aneu <- sweep(ASummary,2,Pmax,FUN="*")

        X <- apply(Pneu,1,range)
        Y <- apply(setValues,1,range)
        colnames(X) <- colnames(Y)
        if (all.equal(X,Y,tolerance=0.01) != TRUE)
            warning("Patterns do not match fixed values.")

        As4fixPs<-list("A"=Aneu,"Asd"=Asd)
        return(As4fixPs)
    }
    else if (setMatrix=="A")
    {
        PSummary <- do.call(cbind,lapply(AP.fixed, function(x) x$Pmean))
        Psd <- do.call(cbind,lapply(AP.fixed, function(x) x$Psd))
        #PSummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Pmean))
        ASummary <- AP.fixed[[1]]$Amean

        Amax <- apply(ASummary,2,max)
        Aneu <- sweep(ASummary,2,Amax,FUN="/")
        Pneu <- sweep(PSummary,1,Amax,FUN="*")

        X <- apply(Aneu,2,range)
        Y <- apply(setValues,2,range)
        colnames(X) <- colnames(Y)
        if (all.equal(X,Y,tolerance=0.01) != TRUE)
            warning("As do not match fixed values.")

        Ps4fixAs<-list("P"=Pneu,"Psd"=Psd)
        return(Ps4fixAs)
    }
    else
    {
        warning("setMatrix can only take values of 'A' or 'P'")
    }

}
