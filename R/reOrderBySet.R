#' reOrderBySet
#'
#' @description <restructures output of gapsRun into a list containing each sets solution for Amean, Pmean, and Asd>
#' @param AP output of gapsRun in parallel
#' @param nFactor number of patterns
#' @param nSets number of sets
#'
#' @return a list containing the \code{nSets} sets solution for Amean under "A", Pmean under "P", and Asd under "Asd"
#' @export
#'
#' @examples \dontrun{
#' reOrderBySet(AP,nFactor,nSets)
#' }
#'
reOrderBySet<-function(AP, #output of gapsRun in parallel
	nFactor, #number of patterns
	nSets #number of sets for parallelization
){
P<-do.call(rbind,lapply(AP, function(x) x$Pmean))
rownames(P)<-paste(rep(1:nSets,each=nFactor),rep(1:nFactor,nSets),sep=".")
A<-lapply(AP, function(x) x$Amean)
Asd<-lapply(AP, function(x) x$Asd)
names(A)=names(Asd)<-paste(rep("Set",nSets),rep(1:nSets),sep="")
return(list("A"=A,"Asd"=Asd,"P"=P))
}
