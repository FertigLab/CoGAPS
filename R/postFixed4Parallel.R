#' postFixed4Parallel
#'
#' @param AP.fixed output of parallel gapsMapRun calls with same FP
#' @param setPs data.frame with rows giving fixed patterns for P used as input for gapsMapRun
#'
#' @return list of two data.frames containing the A matrix estimates or their corresponding standard deviations
#' from output of parallel gapsMapRun
#' @export
#'
#'
postFixed4Parallel <- function(AP.fixed=NA, # output of parallel gapsMapRun calls with same FP
	setPs=NA # data.frame with rows giving fixed patterns for P used as input for gapsMapRun
){

ASummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Amean))
Asd <- do.call(rbind,lapply(AP.fixed, function(x) x$Asd))
#PSummary <- do.call(rbind,lapply(AP.fixed, function(x) x$Pmean))
PSummary <- AP.fixed[[1]]$Pmean

Pmax <- apply(PSummary,1,max)
Pneu <- sweep(PSummary,1,Pmax,FUN="/")
Aneu <- sweep(ASummary,2,Pmax,FUN="*")

X <- apply(Pneu,1,range)
Y <- apply(setPs,1,range)
colnames(X) <- colnames(Y)
if(all.equal(X,Y,tolerance=1e-5)!=TRUE){warning("Patterns do not match fixed values.")}

As4fixPs<-list("A"=Aneu,"Asd"=Asd)
return(As4fixPs)

}
